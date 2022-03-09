
plot_simple_example <- function() {
    library("ggplot2")
    
    ## simulate data
    set.seed(23)
    N <- 100 # number of samples to draw
    mu <- c(-5, 2, 10) # mean of each gaussian
    sd <- c(0.5, 1, 0.25) # standard deviation of each gaussian
    probs <- c(1/3, 1/3, 2/3) # probability of an observation being in each gaussian
    
    df <- data.frame(
        x = sample(x = c(1,2,3), size = N, replace = TRUE, prob = probs),
        y = rnorm(N, mu[x], sd[x])
    )
    
    p <- ggplot()
    
}