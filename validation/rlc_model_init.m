function [model] = rlc_model_init(model, freq_list, range_list)
% generates model of RLC bridge errors for each range and frequency

    % freqs count
    F = numel(freq_list);
    
    % Z axis spots
    N = 100;
    
    % store range list
    model.range_list = range_list;
    
    % extension by virtual zero range
    rngs = [0;range_list]; 
    
    % for each bridge range:
    for r = 1:numel(range_list)
        
        if model.random
            % generate some random gain-phase error
            model.rng(r).f = freq_list;
            model.rng(r).Z = linspace(0,1,100);
            for fid = 1:F
                if strcmpi(model.lin_mode,'1p')
                    % common error for whole range of Z
                    model.rng(r).re_gain(fid,:) = 1 + repmat(linrand(model.gain(1), model.gain(2)), [1 N]);
                    model.rng(r).im_gain(fid,:) = model.rng(r).re_gain(fid,:) + repmat(linrand(model.reim_sym(1), model.reim_sym(2)), [1 N]);
                    model.rng(r).phi(fid,:) = repmat(linrand(model.phase(1), model.phase(2)), [1 N]);
                elseif strcmpi(model.lin_mode,'2p')
                    % linearly changing error between two spots
                    if r == 1
                        gain_a = [0,0];
                        gain_b = [0,0];
                        phi_a = [0,0];
                    else
                        gain_a = model.gain;
                        gain_b = model.reim_sym;
                        phi_a = model.phase;
                    endif                    
                    model.rng(r).re_gain(fid,:) = 1 + interp1([rngs(r) rngs(r+1)], [linrand(gain_a(1),gain_a(2)) linrand(model.gain(1), model.gain(2))], model.rng(r).Z, 'linear', 'extrap');
                    model.rng(r).im_gain(fid,:) = model.rng(r).re_gain(fid,:) + interp1([rngs(r) rngs(r+1)], [linrand(gain_b(1),gain_b(2)) linrand(model.reim_sym(1), model.reim_sym(2))], model.rng(r).Z, 'linear', 'extrap');
                    model.rng(r).phi(fid,:) = interp1([rngs(r) rngs(r+1)], [linrand(phi_a(1),phi_a(1)) linrand(model.phase(1), model.phase(2))], model.rng(r).Z, 'linear', 'extrap');
                elseif strcmpi(model.lin_mode,'poly')
                    % general polynome
                    error('polynomic non-linearity not implemented yet!');                
                endif                 
            endfor                        
            
        else
            % generate fixed gain-phase error
            model.rng(r).f       = freq_list;
            model.rng(r).Z       = linspace(0,1,100);
            model.rng(r).re_gain = repmat(1 + model.gain(1), [F 100]);
            model.rng(r).im_gain = repmat(1 + model.gain(1), [F 100]);
            model.rng(r).phi     = repmat(model.phase(1), [F 100]);
        endif
        
        if ~model.elipse
            % unify re-im gains if no eliptical error enabled 
            model.rng(r).im_gain = model.rng(r).re_gain;                
        endif 
        
    endfor
    
endfunction