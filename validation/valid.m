clear all;
clc;

% this is needed for COM Excel interface
pkg load windows;
pkg load io;

% this folder
mfld = fileparts(mfilename('fullpath'));
cd(mfld);


xls_path = fullfile(mfld,'RLC_fix.xlsx');
xls_copy = fullfile(mfld,'RLC_fix_proc.xlsx');
sheet_refz = 'Ref Z';
sheet_refz_list = 'Ref Z list';
sheet_rng_list = 'Ranges';
sheet_cal = 'Cal Data';
sheet_meas = 'Measurement';


% open source XLS (using COM interface otherwise with OCT it will screw up rest of the sheets!)
xls = xlsopen(xls_path,true,'COM',false);

try % use try-catch to prevent unclosed COM ref

    % get refz sheet dataz
    frz = xls_get_col(xls,sheet_refz,'f');
    Zrn = xls_get_col(xls,sheet_refz,'Nom Z');
    Zrsi = xls_get_col(xls,sheet_refz,'Rs-Xs mult');
    
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
    xls = xls_set_col(xls,sheet_refz,'Rs',Rsr./Zrsi);
    xls = xls_set_col(xls,sheet_refz,'Xs',Xsr./Zrsi);
    xls = xls_set_col(xls,sheet_refz,'U(Rs)',u_Rsr./Zrsi);
    xls = xls_set_col(xls,sheet_refz,'U(Xs)',u_Xsr./Zrsi);
                      
    
    
    % get calibration sheet dataz
    rngc = xls_get_col(xls,sheet_cal,'Range Z',4);
    fcz  = xls_get_col(xls,sheet_cal,'f',4);
    Zcn  = xls_get_col(xls,sheet_cal,'Nom Z',4);     
    Zcsi = xls_get_col(xls,sheet_cal,'Rs-Xs mult',4);
   
    % search refz coresponding to cal. data
    czid = [];
    for k = 1:numel(fcz)
        czid(k,1) = find(fcz(k) == frz & Zcn(k) == Zrn,1);
    endfor
    
    % generate some fake measurements    
    calz_rand_abs = 0.000005;
    calz_rand_rel = 0.000050;
    calz_unc_Rs_abs = 0.000005;
    calz_unc_Rs_min = 0.000005;
    calz_unc_Rs_max = 0.000010;
    calz_unc_Xs_abs = 0.000005;
    calz_unc_Xs_min = 0.000005;
    calz_unc_Xs_max = 0.000010;
    Rsc    = Rsr(czid).*(1 + calz_rand_rel*randn(size(fcz))) + calz_rand_abs*randn(size(fcz));
    ua_Rsc = Rsc.*logrand(calz_unc_Rs_min,calz_unc_Rs_max,size(fcz)) + linrand(0,calz_unc_Rs_abs,size(fcz));
    Xsc    = Xsr(czid).*(1 + calz_rand_rel*randn(size(fcz))) + calz_rand_abs*randn(size(fcz));
    ua_Xsc = Xsc.*logrand(calz_unc_Xs_min,calz_unc_Xs_max,size(fcz)) + linrand(0,calz_unc_Xs_abs,size(fcz));
    
    % calculate calibration factors
    cal_Rs = Rsr(czid) - Rsc;
    cal_Xs = Xsr(czid) - Xsc;
    
    % store calibration data
    xls = xls_set_col(xls,sheet_cal,'Rs',Rsc./Zcsi,4);
    xls = xls_set_col(xls,sheet_cal,'Xs',Xsc./Zcsi,4);
    xls = xls_set_col(xls,sheet_cal,'ua(Rs)',ua_Rsc./Zcsi,4);
    xls = xls_set_col(xls,sheet_cal,'ua(Xs)',ua_Xsc./Zcsi,4);
    
    
    % -- generate measurement data
    
    % unique calibrated ranges
    rngs = unique(rngc);
    % unique calibrated frequencies
    frqs = unique(fcz);
    
    
    % get count of measurements in the sheet (values ignored) 
    fm = xls_get_col(xls, sheet_meas, 'f',4);
    M = numel(fm);
    
    % get measurement scaling factors 
    Zmsi = xls_get_col(xls, sheet_meas, 'Rs-Xs mult',4);
    
    % generate measurement frequencies
    fx   = frqs(round(linrand(1,numel(frqs),[M 1])));
    
    % generate ranges
    rngx = rngs(round(linrand(1,numel(rngs),[M 1])));
    
    % generate test impedances
    %  module |Z| somewhere inside of range
    Zx = [];
    rngs_temp = [0;rngs];
    for k = 1:M
        rid = find(rngx(k) == rngs_temp);
        Zx(k,1) = linrand(max(1e-6,rngs_temp(rid-1)), rngs_temp(rid));
    endfor
    %  random phase angle
    phix = linrand(-pi, +pi, size(Zx));
    %  convert to Rs-Xs
    Rsx = rounddig(Zx.*cos(phix),3);
    Xsx = rounddig(Zx.*sin(phix),3);
    
    % store simulated impedance as a reference
    xls = xls_set_col(xls,sheet_meas,'Valid Rs',Rsx./Zmsi,4);
    xls = xls_set_col(xls,sheet_meas,'Valid Xs',Xsx./Zmsi,4);
    
    
    % distort measurement using inverse calibration data
    Rsm = [];
    Xsm = [];
    for k = 1:M
        % get all matching calibration spots
        cid = find(fx(k) == fcz & rngx(k) == rngc);
        
        if isempty(cid)
            % invalid (this should never happen - wrong range selected or missing calibration data for it)
            
        elseif numel(cid) == 1
            % single spot only (this applies outside boundaries of range)
            
            % apply inverce calibration factor to simulated measurement
            Rsm(k,1) = Rsx(k) - cal_Rs(cid);
            Xsm(k,1) = Xsx(k) - cal_Xs(cid);
            
        else
            % interpolated mode:
            
            % here we have to iteratively adjust |Zx|, because XLS calculates interpolates from measured Z, whereas here we calculate from simulated Z, which is slightly different
            Zx_temp = Zx(k);            
            for r = 1:3    
                % get calibration factor
                c_Rs = interp1(Zcn(cid), cal_Rs(cid), Zx_temp, 'linear', 'extrap');
                c_Xs = interp1(Zcn(cid), cal_Xs(cid), Zx_temp, 'linear', 'extrap');
                
                % apply inverce calibration factor to simulated measurement
                Rsm(k,1) = Rsx(k) - c_Rs;
                Xsm(k,1) = Xsx(k) - c_Xs;
                
                % calculate corrected Zx and recalculate correction
                Zx_temp = (Rsm(k).^2 + Xsm(k).^2).^0.5;                        
            endfor
            
        endif
                
    endfor
    
    % store distorted measurement data
    xls = xls_set_col(xls,sheet_meas,'f',fx,4);
    xls = xls_set_col(xls,sheet_meas,'Range',rngx*1000,4);
    xls = xls_set_col(xls,sheet_meas,'Rs',Rsm./Zmsi,4);
    xls = xls_set_col(xls,sheet_meas,'Xs',Xsm./Zmsi,4);  
    
    
catch err
    try
        xlsclose(xls);
    end               
    error(err.identifier, err.message);
end_try_catch

% always close so there are no zombies of Excel instances!
xlsclose(xls);


return

mcc = 10000;

[xls.num,xls.txt] = xlsread(data_path,data_sheet);
f = csv_get_col(xls,'f',data_row);
w = 2*pi*f;

% capacitor of divider
Cp = csv_get_col(xls,'C',data_row) + randn(1,mcc).*csv_get_col(xls,'UC',data_row);
D  = csv_get_col(xls,'D',data_row) + randn(1,mcc).*csv_get_col(xls,'UD',data_row);
Zc = 1./(w.*Cp.*(j + D));

% resistor of divider
Rs  = csv_get_col(xls,'R',data_row) + randn(1,mcc).*csv_get_col(xls,'UR',data_row);
tau = csv_get_col(xls,'tau',data_row) + randn(1,mcc).*csv_get_col(xls,'Utau',data_row);
Zr = Rs + j*w.*tau.*Rs;

% connection cable and load estimate
Rs = 50e-3 + 50e-3*randn(1,mcc);
Ls = 50e-9 + 50e-9*randn(1,mcc);
Cp = 10e-12 + 20e-12*randn(1,mcc);
Zcs = Rs + j*w*Ls;
Zcp = 1./(j*w.*Cp);

% effective shunt
Zr = 1./(1./Zr + 1./Zcp);

% effective series impedance
Zc = Zc + Zcs;

% tfer
k = Zr./(Zr + Zc);

ratio = mean(abs(k),2)
u_ratio = std(abs(k),[],2)

phase = mean(angle(k),2)
u_phase = std(angle(k),[],2)

kre = mean(real(k),2)
u_kre = std(real(k),[],2)

kim = mean(imag(k),2)
u_kim = std(imag(k),[],2)

% store results
xls = csv_store_col(xls, 'ratio', ratio, data_row);
xls = csv_store_col(xls, 'Uratio', u_ratio, data_row);
xls = csv_store_col(xls, 'phi', 180/pi*phase, data_row);
xls = csv_store_col(xls, 'Uphi', 180/pi*u_phase, data_row);
xls = csv_store_col(xls, 're', kre, data_row);
xls = csv_store_col(xls, 'Ure', u_kre, data_row);
xls = csv_store_col(xls, 'im', kim, data_row);
xls = csv_store_col(xls, 'Uim', u_kim, data_row);

xlswrite(data_path, cat(1,xls.txt,num2cell(xls.num)), [data_sheet '_res']);



