## Simulated data for bayesGDS binary example

if (!require(mvtnorm)) {
    stop("package mvtnorm required")
}
if (!require(devtools)) {
    stop("package devtools required")
}


set.seed(123)

binary_sim <- function(N, k, T) {

    x.mean <- rep(0,k)
    x.cov <- diag(k)
    x.cov[1,1] <- .02
    x.cov[k,k] <- x.cov[1,1]
    mu <- seq(-2,2,length=k)
    Omega <- diag(k)

    X <- t(rmvnorm(N, mean=x.mean, sigma=x.cov)) ## k x N
    B <- t(rmvnorm(N, mean=mu, sigma=Omega)) ## k x N
    XB <- colSums(X * B)
    log.p <- XB - log1p(exp(XB))
    Y <- sapply(log.p, function(q) return(rbinom(1,T,exp(q))))

    binary <- list(Y=Y, X=X, T=T)
    return(binary)
}

binary_small <- binary_sim(20, 2, 200)
devtools::use_data(binary_small, overwrite=TRUE)

binary_large <- binary_sim(800, 3, 300)
devtools::use_data(binary_large, overwrite=TRUE)


