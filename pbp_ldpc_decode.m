function [codeword, success] = pbp_ldpc_decode(H, Pn, MAXITER)
    % H: Parity check matrix
    % Pn: Channel log likelihoods
    % MAXITER: Maximum number of iterations

    [M, N] = size(H);
    
    % Initialization
    P_n_m = zeros(M, N, 2);
    for n = 1:N
        for m = find(H(:, n))' % find the nonzero elements
            % Set P_n_to_m = p_n(b)
            P_n_m(m, n, 1) = Pn(n);
            P_n_m(m, n, 2) = 1 - Pn(n);
        end
    end
    
    for iter = 1:MAXITER
        % Check Node to Variable Node Step
        P_m_n = checkToVariable(H, P_n_m);
        
        % Variable Node to Check Node Step
        [P_n_m, P_n_out] = variableToCheck(H, P_m_n, Pn);
        
        % Make hard decision
        codeword = P_n_out(:, 1) > 0.5;
        
        % Check parity
        if all(mod(H * codeword, 2) == 0)
            success = "Success";
            return;
        end
    end
    
    success = "Fail";
end

function P_m_n = checkToVariable(H, P_n_m)
    [M, N, ~] = size(P_n_m);
    P_m_n = zeros(M, N, 2);
    
    for m = 1:M
        for n = find(H(m, :))
            neighbors = find(H(m, :));
            neighbors(neighbors == n) = [];
            
            delta = prod(P_n_m(m, neighbors, 2) - P_n_m(m, neighbors, 1));
            P_m_n(m, n, 1) = (1 - delta) / 2;
            P_m_n(m, n, 2) = (1 + delta) / 2;
        end
    end
end

function [P_n_m, P_n_out] = variableToCheck(H, P_m_n, Pn)
    [M, N, ~] = size(P_m_n);
    P_n_m = zeros(M, N, 2);
    P_n_out = zeros(N, 2);
    
    for n = 1:N
        for m = find(H(:, n))'
            neighbors = find(H(:, n));
            neighbors(neighbors == m) = [];
            
            P_n_m(m, n, 1) = Pn(n) * prod(P_m_n(neighbors, n, 1));
            P_n_m(m, n, 2) = (1 - Pn(n)) * prod(P_m_n(neighbors, n, 2));
            
            alpha = 1 / (P_n_m(m, n, 1) + P_n_m(m, n, 2));
            P_n_m(m, n, :) = alpha * P_n_m(m, n, :);
        end
        out_neighbors = find(H(:, n));
        P_n_out(n, 1) = Pn(n) * prod(P_m_n(out_neighbors, n, 1));
        P_n_out(n, 2) = (1 - Pn(n)) * prod(P_m_n(out_neighbors, n, 2));
        
        alpha = 1 / (P_n_out(n, 1) + P_n_out(n, 2));
        P_n_out(n, :) = alpha * P_n_out(n, :);
    end
end