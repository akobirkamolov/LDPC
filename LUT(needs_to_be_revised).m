function [T] = LUT(M, N, delta)
% Initialize the look-up table
T = zeros(M, N);

% Calculate the values at each index
for r = 1:M
    for c = 1:N
        % Box-plus inputs
        a = (r-1)*delta + delta/2;
        b = (c-1)*delta + delta/2;

        % Box-plus
        check = 2*atanh(tanh(a/2)*tanh(b/2));

        % Quantize and assign
        T(r,c) = abs(floor(check/delta));
    end
end
end