# variational inference --------------------------------------------------------
vi_dirichlet_mixture <- function(
    y,
    K = 20,
    priors = NULL,
    max_iter = 30
) {
    
    if (is.null(priors)) {
        lambda = 0.1
        mu_mu_k = 0
        sigmasq_mu_k = 100
        alpha_sigmasq_k = 0.01
        beta_sigmasq_k = 0.01
        a_lambda = 0.01
        b_lambda = 0.01
    } else {
        
        if (!is.list(priors)) {
            stop("priors must be a list")
        }
    }
    
    N <- length(y) # sample size
    sum_y <- sum(y)
    elbos <- rep(-Inf, max_iter)
    
    ## initialize parameters
    tilde_p_q_p <- matrix(0, nrow = N, ncol = K)
    p_q_p <- matrix(0, nrow = N, ncol = K)
    
    mu_q_mu <- rnorm(K, mean = mean(y), sd = sd(y))
    sigmasq_q_mu <- runif(K)
    
    A_q_sigmasq <- runif(K)
    B_q_sigmasq <- runif(K)
    
    alpha_q_V <- runif(K-1)
    beta_q_V <- runif(K-1)
    
    for(m in 1:max_iter) {
        
        ## expectation of sigma^2
        E_sigma_sq <- B_q_sigmasq/A_q_sigmasq
        
        ## optimal density for z_ik
        for(k in 1:K) {
            ### stick-breaking representation
            E_p_k <- 0
            if(k == 1) {
                E_p_k <- digamma(alpha_q_V[k]) - digamma(alpha_q_V[k] + beta_q_V[k])
            } else if (k > 1 & k < K) {
                E_p_k <- digamma(alpha_q_V[k]) - digamma(alpha_q_V[k] + beta_q_V[k]) + 
                    sum(digamma(beta_q_V[1:(k-1)]) - 
                            digamma(alpha_q_V[1:(k-1)] + beta_q_V[1:(k-1)]))
            } else { # if k == K then E[log V_K] = E[log 1] = 0
                E_p_k <- sum(digamma(beta_q_V[1:(k-1)]) - 
                                  digamma(alpha_q_V[1:(k-1)] + beta_q_V[1:(k-1)]))
            }
            
            tilde_p_q_p[, k] <- E_p_k - 
                (1/2)*(digamma(A_q_sigmasq[k]) + log(B_q_sigmasq[k]) + 
                           (1/E_sigma_sq)*((y-mu_q_mu[k])^2 + sigmasq_q_mu[k]))
        }
        
        p_q_p <- t(apply(tilde_p_q_p, 1, function(x)exp(x)/sum(exp(x))))
        
        E_n <- colSums(p_q_p) # a count of observations in each group
        
        for(k in 1:K){
            
            ## optimal density for V_k
            if(k < K) {
                alpha_q_V[k] <- 1 + E_n[k]
                beta_q_V[k] <- lambda + sum(E_n[(k+1):K])
            }
            
            ## optimal density for mu_k
            sigmasq_q_mu[k] <- 1/(1/sigmasq_mu_k + (1/E_sigma_sq[k])*E_n[k])
            mu_q_mu[k] <- sigmasq_q_mu[k]*(mu_mu_k/sigmasq_mu_k + (1/E_sigma_sq[k])*sum(y*p_q_p[, k]))
            
            ## optimal density for sigma^2_k
            A_q_sigmasq[k] <- alpha_sigmasq_k + E_n[k]/2
            B_q_sigmasq[k] <- beta_sigmasq_k + (1/2)*(sum(((y - mu_q_mu[k])^2 + sigmasq_q_mu[k])*p_q_p[, k]))
        }
        
    }
    
    return(
        list(
            p_q_p = p_q_p,
            alpha_q_V = alpha_q_V,
            beta_q_V = beta_q_V,
            mu_q_mu = mu_q_mu,
            sigmasq_q_mu = sigmasq_q_mu,
            A_q_sigmasq = A_q_sigmasq,
            B_q_sigmasq = B_q_sigmasq,
            elbos = elbos
        )
    )
}