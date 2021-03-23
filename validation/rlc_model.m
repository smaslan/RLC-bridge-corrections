function [Rs_out, Xs_out] = rlc_model(model, f, range, Rs, Xs, Zmod, inverse)
% applies previously initialized RLC bridge model to measured data

    % impedance modulus
    Z = complex(Rs, Xs);
    
    % reference |Z| for interpolation
    if isempty(Zmod)
        Zmod = abs(Z);
    else
        Zmod = abs(Zmod);
    endif
    
    % for each input:
    err_gain = [];
    err_phi = [];
    for r = 1:numel(range)
        
        % get range from model ranges list
        rid = find(range(r) == model.range_list, 1);
        if isempty(rid)
            error(sprintf('Range %g not found in model!',range(r)));
        endif
        rng = model.rng(rid);
        %rng = model.rng(1);
        
        % interpolate RLC errors model to get error at this spot
        err_re_gain(r,1) = interp2(rng.Z, rng.f, rng.re_gain, Zmod(r), f(r), 'linear');
        err_im_gain(r,1) = interp2(rng.Z, rng.f, rng.im_gain, Zmod(r), f(r), 'linear');
        err_phi(r,1)  = interp2(rng.Z, rng.f, rng.phi,  Zmod(r), f(r), 'linear');
                
    endfor
    
    % apply gain-phase
    if inverse
        Z_out = complex(real(Z)./err_re_gain, imag(Z)./err_im_gain).*exp(-j*err_phi);
    else                
        Z_out = complex(real(Z).*err_re_gain, imag(Z).*err_im_gain).*exp(+j*err_phi);
    endif
    
    % return components
    Rs_out = real(Z_out);
    Xs_out = imag(Z_out);
    
endfunction