function []=writeout(varargin)
% -writeout(fid,array,header) 
% quickly write out a data array with a header;
%- fid, file's handle
% -array, data array
% -header,strings

if nargin == 2
    fid = varargin{1};
    data = varargin{2};
    
    [m,n]=size(data);
    for i=1:m
        for j=1:n
            if j==n
                fprintf(fid,'%10.6f\n',data(i,j));
            else
                fprintf(fid,'%10.6f\t',data(i,j));
            end
        end
    end
    fclose(fid);
    
elseif nargin == 3
    fid = varargin{1};
    data = varargin{2};
    header = varargin{3}; % the inputed header has to be a string
    fprintf(fid,'%s\n',header);
    
    [m,n]=size(data);
    for i=1:m
        for j=1:n
            if j==n
                fprintf(fid,'%10.6f\n',data(i,j));
            else
                fprintf(fid,'%10.6f\t',data(i,j));
            end
        end
    end
    fclose(fid);
    
    
elseif nargin == 4
    
    fid = varargin{1};
    data = varargin{2};
    header = varargin{3}; % the inputed header has to be a string
    fprintf(fid,'%s\n',header);
    ID = varargin{4};
    
    [m,n]=size(data);
    for i=1:m
        for j=1:n
            if j==n
                fprintf(fid,'%10.6f  %s\n',data(i,j),char(ID(i)));
            else
                fprintf(fid,'%10.6f\t',data(i,j));
            end
        end
    end
    fclose(fid);
    
else
    error('wrong number of input variant')
end

end

