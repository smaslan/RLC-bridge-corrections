%====================================================================================================
% This is simple script that validates the RLC bridge corrections XLS spreadsheet file.
% It generates simulated reference standards, simulated calibration measurements of RLC and 
% then simulates measurements with errors of RLC meter. If XLS works right, the XLS sheet should
% shown no deviation from simulated data.
%
% (c) 2021 Stanislav Maslan, s.maslan@seznam.cz.
% The script is distributed under MIT license, https://opensource.org/licenses/MIT. 
%====================================================================================================
clear all;
clc;

% this is needed for COM Excel interface
pkg load windows;
pkg load io;

% this folder
mfld = fileparts(mfilename('fullpath'));
cd(mfld);

% correction XLS to validate
xls_path = fullfile(mfld,'../spreadsheets/RLC_fix.xlsx');
% sheet name with reference impedances
sheet_refz = 'Ref Z';
% sheet name with list of used reference impedances
sheet_refz_list = 'Ref Z list';
% sheet name with calibration data
sheet_cal = 'Cal Data';
% sheet name to which it will simulate measurements
sheet_meas = 'Measurement';

% open source XLS (using COM interface otherwise with OCT it will screw up rest of the sheets!)
xls = xlsopen(xls_path,true,'COM',false);

try % usin try-catch to prevent unclosed XLS COM ref

    % -- generate reference Z data:
    % note: XLS must contain tables with all frequencies and ranges and nominal Z values. This code will only fill in the sheet.
    
    % get refz sheet dataz
    frz = xls_get_col(xls,sheet_refz,'f',4);
    Zrn = xls_get_col(xls,sheet_refz,'Nom Z',4);
    Zrsi = xls_get_col(xls,sheet_refz,'Rs-Xs mult',4);
    
    % generate some reference impedance values 
    refz_rand = 0.0001;
    refz_fdep_Rs_max_dev = linrand(0.00,0.1);
    refz_fdep_Rs_pow = 1.5;
    refz_unc_Rs_max = 0.000500;
    refz_unc_Rs_min = 0.000010;
    refz_unc_Rs_pow = 1.5;
    refz_unc_Xs_max = 0.000500;
    refz_unc_Xs_min = 0.000050;
    refz_unc_Xs_pow = 1.0;
    refz_fdep_Rs_pow = 1.5;
    refz_fdep_Xs_max_dev = linrand(-0.05,0.05);
    refz_fdep_Xs_pow = 1.0;
    Rsr   = Zrn.*(1 + refz_rand*randn(size(Zrn)) + (frz/max(frz)).^refz_fdep_Rs_pow*refz_fdep_Rs_max_dev);
    u_Rsr = Zrn.*(max((frz/max(frz)).^refz_unc_Rs_pow*refz_unc_Rs_max,refz_unc_Rs_min));
    Xsr   = Zrn.*(0 + refz_rand*randn(size(Zrn)) + (frz/max(frz)).^refz_fdep_Xs_pow*refz_fdep_Xs_max_dev);
    u_Xsr = Zrn.*(max((frz/max(frz)).^refz_unc_Xs_pow*refz_unc_Xs_max,refz_unc_Xs_min));
    
    % store Ref Z data
    xls = xls_set_col(xls,sheet_refz,'Rs',Rsr./Zrsi,4);
    xls = xls_set_col(xls,sheet_refz,'Xs',Xsr./Zrsi,4);
    xls = xls_set_col(xls,sheet_refz,'U(Rs)',u_Rsr./Zrsi,4);
    xls = xls_set_col(xls,sheet_refz,'U(Xs)',u_Xsr./Zrsi,4);



    % -- generate calibration data:
    % note: XLS must contain tables with all frequencies and ranges and nominal Z values. This code will only fill in the sheet.
    
    
    % Definining bridge errors model:
    model = struct();
    %   mode of error (gain-phase)
    model.lin_mode = '2p'; 
    %   min max gain error array (deviation from 1.000000)
    model.gain   = [-0.000100 +0.000100];
    %   enable different gain for real-imag axes? 
    model.elipse = 0;
    %   elipticity of complex plane 
    model.reim_sym = [-0.000010 +0.000010];
    %   min max phase error [rad]
    model.phase  = [-0.000100 +0.000100];
    %   randomize? if 0, will generate constant errors 
    model.random = 1; 
                              
   
    % get calibration sheet dataz
    rngc = xls_get_col(xls,sheet_cal,'Range Z',4);
    fcz  = xls_get_col(xls,sheet_cal,'f',4);
    Zcn  = xls_get_col(xls,sheet_cal,'Nom Z',4);     
    Zcsi = xls_get_col(xls,sheet_cal,'Rs-Xs mult',4);
    
    % get all existing calibration standads nominals
    std_list = unique(Zcn);
    
    % get all existing ranges
    rng_list = unique(rngc);
    
    % get all existing calibration frequencies
    freq_list = unique(fcz);
    
    % initialize RLC bridge errors model
    model = rlc_model_init(model, freq_list, rng_list);
   
    % search refz coresponding to cal. data
    czid = [];
    for k = 1:numel(fcz)
        czid(k,1) = find(fcz(k) == frz & Zcn(k) == Zrn,1);
    endfor
       
    % generate some fake measurements    
    calz_unc_abs = 0.000002;
    calz_unc_min = 0.000005;
    calz_unc_max = 0.000010;
    % apply RLC bridge model
    [Rsc, Xsc] = rlc_model(model, fcz, rngc, Rsr(czid), Xsr(czid), [], 0);
    % generate some uncertainties
    ua_Rsc = Rsc.*logrand(calz_unc_min,calz_unc_max,size(fcz)) + linrand(0,calz_unc_abs,size(fcz));    
    ua_Xsc = Xsc.*logrand(calz_unc_min,calz_unc_max,size(fcz)) + linrand(0,calz_unc_abs,size(fcz));
    
    % generate short resodual for cal sheet
    calz_sh_rnd = 0.000002;
    Rsc_sh = linrand(-calz_sh_rnd, calz_sh_rnd, size(fcz));
    Xsc_sh = linrand(-calz_sh_rnd, calz_sh_rnd, size(fcz));
    ua_Rsc_sh = Rsc_sh.*logrand(calz_unc_min,calz_unc_max,size(fcz)) + linrand(0,calz_unc_abs,size(fcz));    
    ua_Xsc_sh = Xsc_sh.*logrand(calz_unc_min,calz_unc_max,size(fcz)) + linrand(0,calz_unc_abs,size(fcz));
    Rsc_sh(Zcn == 0) = NaN;
    Xsc_sh(Zcn == 0) = NaN;
    ua_Rsc_sh(Zcn == 0) = NaN;
    ua_Xsc_sh(Zcn == 0) = NaN;
        
    % store calibration data
    xls = xls_set_col(xls,sheet_cal,'Rs',(Rsc + Rsc_sh)./Zcsi,4);
    xls = xls_set_col(xls,sheet_cal,'Xs',(Xsc + Xsc_sh)./Zcsi,4);
    xls = xls_set_col(xls,sheet_cal,'ua(Rs)',ua_Rsc./Zcsi,4);
    xls = xls_set_col(xls,sheet_cal,'ua(Xs)',ua_Xsc./Zcsi,4);
    xls = xls_set_col(xls,sheet_cal,'sh Rs',Rsc_sh./Zcsi,4);
    xls = xls_set_col(xls,sheet_cal,'sh Xs',Xsc_sh./Zcsi,4);
    xls = xls_set_col(xls,sheet_cal,'sh ua(Rs)',ua_Rsc_sh./Zcsi,4);
    xls = xls_set_col(xls,sheet_cal,'sh ua(Xs)',ua_Xsc_sh./Zcsi,4);
    
    
    
    % -- generate measurement data:
    % note: XLS must contain some rows to fill in. It will simulate as much measurements as valid rows find
    
    % get count of measurements in the sheet (values ignored) 
    fm = xls_get_col(xls, sheet_meas, 'f',4);
    M = numel(fm);
    
    % get measurement scaling factors 
    Zmsi = xls_get_col(xls, sheet_meas, 'Rs-Xs mult',4);
    
    % generate measurement frequencies
    fx   = freq_list(round(linrand(1,numel(freq_list),[M 1])));
    
    % generate ranges
    rngx = rng_list(round(linrand(1,numel(rng_list),[M 1])));
    
    % generate test impedances
    %  over-ranging (0.05: up to +-5% outside optimal range) 
    measz_over = 0.05;
    %  module |Z| somewhere inside of range
    Zx = [];
    rngs_temp = [0;rng_list];
    for k = 1:M
        rid = find(rngx(k) == rngs_temp);
        delta = rngs_temp(rid) - rngs_temp(rid-1);        
        Zx(k,1) = linrand(max(1e-6,rngs_temp(rid-1) - delta*measz_over), rngs_temp(rid) + delta*measz_over);
    endfor
    %  random phase angle
    phix = linrand(-pi, +pi, size(Zx));
    %  convert to Rs-Xs
    Rsx = rounddig(Zx.*cos(phix),3);
    Xsx = rounddig(Zx.*sin(phix),3);
    
    % generate short residual for measurement sheet
    measz_sh_rnd  = 0.000002;
    measz_unc_abs = 0.000002;
    measz_unc_min = 0.000005;
    measz_unc_max = 0.000010;
    Rsx_sh = linrand(-measz_sh_rnd, measz_sh_rnd, size(fx));
    Xsx_sh = linrand(-measz_sh_rnd, measz_sh_rnd, size(fx));
    ua_Rsx_sh = Rsx_sh.*logrand(measz_unc_min,measz_unc_max,size(fx)) + linrand(0,measz_unc_abs,size(fx));    
    ua_Xsx_sh = Xsx_sh.*logrand(measz_unc_min,measz_unc_max,size(fx)) + linrand(0,measz_unc_abs,size(fx));    
    ua_Rsx = Rsx.*logrand(measz_unc_min,measz_unc_max,size(fx)) + linrand(0,measz_unc_abs,size(fx));    
    ua_Xsx = Xsx.*logrand(measz_unc_min,measz_unc_max,size(fx)) + linrand(0,measz_unc_abs,size(fx));
       
    
    % store simulated impedance as a reference
    xls = xls_set_col(xls,sheet_meas,'Valid Rs',Rsx./Zmsi,4);
    xls = xls_set_col(xls,sheet_meas,'Valid Xs',Xsx./Zmsi,4);
    
    % apply RLC bridge error to the data (iterative fix because XLS uses measured |Z| for interpolation, and we have actual |Z| only)
    Rsm = Rsx; Xsm = Xsx;
    for it = 1:5
        Ztmp = (Rsm.^2 + Xsm.^2).^0.5;
        [Rsm, Xsm] = rlc_model(model, fx, rngx, Rsx, Xsx, Ztmp, 0);        
    endfor
         
    % store distorted measurement data
    xls = xls_set_col(xls,sheet_meas,'f',fx,4);
    xls = xls_set_col(xls,sheet_meas,'Range',rngx*1000,4);
    xls = xls_set_col(xls,sheet_meas,'Rs',(Rsm + Rsx_sh)./Zmsi,4);
    xls = xls_set_col(xls,sheet_meas,'Xs',(Xsm + Xsx_sh)./Zmsi,4);  
    xls = xls_set_col(xls,sheet_meas,'ua(Rs)',(ua_Rsx)./Zmsi,4);
    xls = xls_set_col(xls,sheet_meas,'ua(Xs)',(ua_Xsx)./Zmsi,4);
    xls = xls_set_col(xls,sheet_meas,'sh Rs',(Rsx_sh)./Zmsi,4);
    xls = xls_set_col(xls,sheet_meas,'sh Xs',(Xsx_sh)./Zmsi,4);
    xls = xls_set_col(xls,sheet_meas,'sh ua(Rs)',(ua_Rsx_sh)./Zmsi,4);
    xls = xls_set_col(xls,sheet_meas,'sh ua(Xs)',(ua_Xsx_sh)./Zmsi,4);
    
    
catch err
    try
        xlsclose(xls);
    end               
    error(err);
end_try_catch

% always close so there are no zombies of Excel instances!
xlsclose(xls);



