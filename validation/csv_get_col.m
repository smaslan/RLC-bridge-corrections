function [data,col] = csv_get_col(csv,name,first_data_row=3)
    if isstruct(csv)
        % try to find column by name
        col = find(strcmpi(csv.txt(1,:),name),1);
        if isempty(col)
            error(sprintf('Cannot find column ''%s'' in given XLS data!',name));
        endif
        data = [csv.num(first_data_row:end,col)](:);
    else
        % try to find column by name
        col = find(strcmpi(csv(1,:),name),1);
        if isempty(col)
            error(sprintf('Cannot find column ''%s'' in given CSV data!',name));
        endif
        data = [csv{first_data_row:end,col}](:);
    endif
endfunction