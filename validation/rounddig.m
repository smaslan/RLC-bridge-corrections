function y = rounddig(x,d)    
    digits = ceil(log10(abs(x)));    
    round_base = 10.^-(digits - d);    
    y = round(x.*round_base)./round_base;
    y(x == 0) = 0;
end