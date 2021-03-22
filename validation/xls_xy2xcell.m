function [str]=xls_xy2xcell(x,y)
    
    chr = '';
    if x > 26 && x <= (26+26)
        x -= 26;
        chr = 'A'; 
    elseif x > (26+26)
        x -= 26;
    endif 
    str = sprintf('%s%s%d',chr,dec2base(x-1,'ABCDEFGHIJKLMNOPQRSTUVWXYZ'),y);    
    
endfunction