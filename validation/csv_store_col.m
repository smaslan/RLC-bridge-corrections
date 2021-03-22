function [csv] = csv_store_col(csv,name,data,first_data_row=3)
    if isstruct(csv)
        % try to find column by name
        col = find(strcmpi(csv.txt(1,:),name),1);
        if isempty(col)
            error(sprintf('Cannot find column ''%s'' in given XLS data!',name));
        endif        
        
        csv.num(first_data_row:first_data_row - 1 + numel(data),col) = data;
                
    else
        error('not implemented for raw CSV');
    endif
endfunction