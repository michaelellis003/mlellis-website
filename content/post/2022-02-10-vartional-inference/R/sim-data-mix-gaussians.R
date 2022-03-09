## simulate data
set.seed(23)
N <- 100 # number of samples to draw
mu <- c(-5, 2, 10) # mean of each gaussian
sd <- c(0.5, 1, 0.25) # standard deviation of each gaussian
probs <- c(1/3, 1/3, 2/3) # probability of an observation being in each gaussian

x <- sample(x = c(1,2,3), size = N, replace = TRUE, prob = probs)
y <- rnorm(N, mu[x], sd[x])

source("R/vi_dirichlet_mixture.R")
output <- vi_dirichlet_mixture(
    y,
    K = 20,
    priors = NULL,
    max_iter = 100
) 

round(output$mu_q_mu, 2)
