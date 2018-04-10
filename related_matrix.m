function [S] = related_matrix(code_lenght,related_number)
S = zeros(code_lenght,code_lenght);
for i = 1:code_lenght
    for j = 1:code_lenght
        if i<j
            S(i,j) = related_number^(j-i);
        else
            S(i,j) = related_number^(i-j);
        end
    end
end
end
    