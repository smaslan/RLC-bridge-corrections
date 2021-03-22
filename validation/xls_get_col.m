function [data,col] = xls_get_col(xls,sheet,name,header_row=3,data_row_offset=2)

    % get sheet data
    xdata = xls2oct(xls,sheet);
    
    % try to find column by name
    col = find(strcmpi(xdata(header_row,:),name),1);    
    if isempty(col)
        error(sprintf('Cannot find column ''%s'' in given XLS data!',name));
    endif    
    data = [xdata{header_row+data_row_offset:end,col}](:);

endfunction