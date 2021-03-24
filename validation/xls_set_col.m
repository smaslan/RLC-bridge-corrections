function [xls, col] = xls_set_col(xls,sheet,name,data,header_row,data_row_offset=2)
    
    % get sheet data
    xdata = xls2oct(xls,sheet);
    
    if ischar(header_row)
        % try to search table header row  
        hid = find(strcmpi(xdata(:),header_row));
        if isempty(hid)
            error(sprintf('Octave header tag ''%s'' not found in the XLS sheet ''%s''!',header_row,sheet));
        endif
        header_row = mod(hid,size(xdata,1));
    endif 

    % try to find column by name
    col = find(strcmpi(xdata(header_row,:),name),1);
    if isempty(col)
        error(sprintf('Cannot find column ''%s'' in given XLS data!',name));
    endif
    % write column
    xls = oct2xls(num2cell(data),xls,sheet,[xls_xy2xcell(col,header_row+data_row_offset) ':' xls_xy2xcell(col,header_row+data_row_offset-1+size(data,1))]);
    
endfunction