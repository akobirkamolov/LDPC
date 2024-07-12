function r_hat = BitFlipDecoder(r, H, max_iterations, threshold)
    % r: received vector
    % H: parity check matrix
    % max_iterations: maximum number of iterations
    % threshold: threshold for bit flipping

    % Step 1: Initialize r_hat
    r_hat = r;
    
    for iteration = 1:max_iterations
        % Step 2: Compute syndrome
        s = mod(r_hat * H', 2);
        
        % If syndrome is zero, stop
        if all(s == 0)
            return;
        end
        
        % Step 3: Compute f
        f = sum(H(s == 1, :), 1);
        
        % Step 4: Flip bits
        bits_to_flip = f > threshold;
        r_hat(bits_to_flip) = 1 - r_hat(bits_to_flip);
    end
end