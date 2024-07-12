function [T] = LUT_simple(N)
T = zeros(N);
for i= 1:N
    for j = 1:N
        if(i==j && i == 1)
            T(i,j) = 0;
        elseif(i==j)
            T(i,j) = i-2;
        elseif(i>j)
            T(i,j) = chk_quantizer(i-1,j-1);
        end
    end
end
for i=1:N
    for j=1:N
        if(j>i)
            T(i,j) = T(j,i);
        end
    end
end
end