function [c_hat, status] = bp_ldpc_decode(H, Ln, max_iter)
    [M, N] = size(H);

    % Initialzie messages 
    L_nm = repmat(Ln, M, 1);
    L_mn = zeros(M, N);

    for iter = 1:max_iter
        % Check Node to Variable Node Step
        for m = 1:M
            n_indices = find(H(m, :));
            for n = n_indices
                product = prod(tanh(L_nm(m, H(m, :) == 1) / 2));
                L_mn(m, n) = 2 * atanh(product / tanh(L_nm(m, n) / 2));
            end
        end
        for n = 1:N
            m_indices = find(H(:, n));
            for m = m_indices'
                L_nm(m, n) = Ln(n) + sum(L_mn(H(:, n) == 1, n)) - L_mn(m, n);
            end 
        end

        % Compute output likelihoods 
        Ln_out = Ln + sum(L_mn, 1);

        % Make decisions 
        c_hat = (Ln_out < 0);

        % Check parity
        if all(mod(H * c_hat', 2) == 0)
            status = 'Success';
            return;
        end
    end

    % If max iterations reached wihtout convergence 
    status = 'Fail';
end