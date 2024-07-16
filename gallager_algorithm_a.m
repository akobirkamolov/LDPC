function decoded_bits = gallager_algorithm_a(y, H)
    % Main function to implement Gallager's Algorithm A
    % received_bits: The channel input yn
    % H: Parity check matrix

    B_n_to_m = initialize_B(y, H);
    
    max_iterations = 1;  % Set a maximum number of iterations
    for iter = 1:max_iterations
        B_m_to_n = check_to_variable_node(B_n_to_m, H);
        B_n_to_m = variable_node_computation(B_n_to_m, B_m_to_n, y, H);
    end
    
    decoded_bits = decoding_decision(B_n_to_m, H);
end

function B_n_to_m = initialize_B(y, H)
    % Initialization step
    [M, N] = size(H);
    B_n_to_m = zeros(M, N);
    for n = 1:N
        B_n_to_m(:, n) = y(n) * H(:, n);
    end
end

function B_m_to_n = check_to_variable_node(B_n_to_m, H)
    % Check Node to Variable Node step
    M = size(H, 1);
    B_m_to_n = zeros(size(H));
    for m = 1:M
        neighbors = find(H(m, :));
        for n = neighbors
            B_m_to_n(m, n) = mod((sum(B_n_to_m(m, neighbors))-B_n_to_m(m,n)), 2);
        end
    end
end

function B_n_to_m = variable_node_computation(B_n_to_m, B_m_to_n, y, H)
    % Variable Node Computation step
    [M, N] = size(H);
    for n = 1:N
        for m = 1:M
            neighbors = find(H(:, n));
            neighbors(neighbors == m) = []; 
            if all(B_m_to_n(neighbors, n) ~= y(n))
                if H(m,n)
                    B_n_to_m(m, n) = 1- y(n);
                end
            else
                if H(m,n)
                    B_n_to_m(m, n) = y(n);
                end
            end
        end
    end
end

function decoded_bits = decoding_decision(B_n_to_m, H)
    % Decoding Decision step
    N = size(H, 2);
    decoded_bits = zeros(1, N);
    for n = 1:N
        neighbors = find(H(:, n));
        if mod(length(neighbors), 2) == 0
            decoded_bits(n) = mode([B_n_to_m(neighbors, n); yn(n)]);
        else
            decoded_bits(n) = mode(B_n_to_m(neighbors, n));
        end
    end
end