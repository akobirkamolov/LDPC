import numpy as np

def ldpc_decode(H, Ln, max_iter):
    N = H.shape[1]  # number of variable nodes
    M = H.shape[0]  # number of check nodes
    
    # Initialize messages
    L_nm = np.tile(Ln, (M, 1))
    L_mn = np.zeros((M, N))
    
    for iter in range(max_iter):
        # Check Node to Variable Node Step
        for m in range(M):
            for n in np.where(H[m, :])[0]:
                prod = np.prod(np.tanh(L_nm[m, H[m, :] == 1] / 2))
                L_mn[m, n] = 2 * np.arctanh(prod / np.tanh(L_nm[m, n] / 2))
        
        # Variable Node to Check Node Step
        for n in range(N):
            for m in np.where(H[:, n])[0]:
                L_nm[m, n] = Ln[n] + np.sum(L_mn[H[:, n] == 1, n]) - L_mn[m, n]
        
        # Compute output likelihoods
        Ln_out = Ln + np.sum(L_mn[:, :], axis=0)
        
        # Make decisions
        c_hat = (Ln_out < 0).astype(int)
        
        # Check parity
        if np.all(np.dot(H, c_hat) % 2 == 0):
            return c_hat, "Success"
    
    # If max iterations reached without convergence
    return c_hat, "Fail"