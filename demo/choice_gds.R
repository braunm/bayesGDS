## Runs the hierarchical binary choice example

library(Matrix)
library(trustOptim)
library(sparseHessianFD)
library(sparseMVN)
library(plyr)
library(mvtnorm)

library(doParallel)
run.par <- FALSE
if(run.par) registerDoParallel(cores=2) else registerDoParallel(cores=1)

seed.id <- 123
set.seed(seed.id)


rmvn.sparse.wrap <- function(n.draws, params) {
## sample MVN with sparse precision matrix

    res <- rmvn.sparse(n.draws, params$mean, params$CH, prec=TRUE)
    return(res)
}

dmvn.sparse.wrap <- function(d, params) {
## MVN density with sparse precision

    res <- dmvn.sparse(d, params$mean, params$CH, prec=TRUE)
    return(res)
}


## parameters for simulating data
N <- 100  ## number of heterogeneous units
k <- 2  ## number of covariates
T <- 40  ## number of "purchase opportunities per unit
nvars <- N*k + k

## Simulate data and set priors

x.mean <- rep(0,k)
x.cov <- diag(k)
mu <- rnorm(k,0,10)
Omega <- diag(k)
inv.Sigma <- rWishart(1,k+5,diag(k))[,,1]
inv.Omega <- solve(Omega)
X <- t(rmvnorm(N, mean=x.mean, sigma=x.cov)) ## k x N
B <- t(rmvnorm(N, mean=mu, sigma=Omega)) ## k x N
XB <- colSums(X * B)
log.p <- XB - log1p(exp(XB))
Y <- apply(as.matrix(log.p), 1,function(q) return(rbinom(1,T,exp(q))))

start <- rnorm(nvars) ## random starting values
hess.struct <- bayesGDS::demo.get.hess.struct(N, k)

## Setting up function to compute Hessian using sparseHessianFD package.
obj <- new.sparse.hessian.obj(start, fn=bayesGDS::demo.get.f,
                              gr=bayesGDS::demo.get.grad,
                              hs=hess.struct, Y=Y, X=X,
                              inv.Omega=inv.Omega,
                              inv.Sigma=inv.Sigma, T=T)

get.f.wrap <- function(x)  return(obj$fn(x))
get.df.wrap <- function(x)  return(obj$gr(x))
get.hessian.wrap <- function(x)  return(obj$hessian(x))

cat("Using Sparse method in trust.optim\n")
cat("to find posterior mode\n")

opt <- trust.optim(start, fn=get.f.wrap,
                   gr = get.df.wrap,
                   hs = get.hessian.wrap,
                   method = "Sparse",
                   control = list(
                     start.trust.radius=5,
                     stop.trust.radius = 1e-7,
                     prec=1e-7,
                     function.scale.factor=-1,
                     report.freq=1L,
                     report.level=4L,
                     report.precision=1L,
                     maxit=500L,
                     preconditioner=0L
                     ) 
                   )

post.mode <- opt$solution
hess <- opt$hessian
var.names <- names(post.mode)

## parameters for the GDS algorithm itself
n.draws <- 20  ## total number of draws needed
batch.size <- 10  ## draws per instance of sample.GDS
M <- 20000  ## proposal draws
max.tries <- 100000  ## to keep sample.GDS from running forever
ds.scale <- 0.97  ## scaling factor for proposal density
fn.dens.prop <- dmvn.sparse.wrap
fn.draw.prop<- rmvn.sparse.wrap


chol.hess <- Cholesky(-ds.scale*hess)

prop.params <- list(mean = post.mode,
                    CH = chol.hess
                    )

log.c1 <- opt$fval
log.c2 <- dmvn.sparse.wrap(post.mode, prop.params)

cat("Collecting GDS Proposal Draws\n")
draws.m <- as(fn.draw.prop(M,prop.params),"matrix")
log.post.m <- aaply(draws.m, 1,get.f.wrap, .parallel=run.par)
log.prop.m <- fn.dens.prop(draws.m,params=prop.params)
log.phi <- log.post.m - log.prop.m +log.c2 - log.c1

invalid.scale <- any(log.phi>0)
cat("Are any log.phi > 0?  ",invalid.scale,"\n")

## if invalid.scale is TRUE, need to change
## the proposal density


if (!invalid.scale) {
    n.batch <- floor(n.draws / batch.size)
    
    cat("Generating DS draws - accept-reject phase\n")
    draws.list <- foreach(i=1:n.batch, .inorder=FALSE) %dopar% sample.GDS(
                                      n.draws=batch.size,
                                      log.phi=log.phi,
                                      post.mode=post.mode,
                                      fn.dens.post = get.f.wrap,
                                      fn.dens.prop = dmvn.sparse.wrap,
                                      fn.draw.prop = rmvn.sparse.wrap,
                                      prop.params = prop.params,
                                      max.tries=max.tries,
                                      report.freq=50,
                                      announce=TRUE,
                                      thread.id=i,
                                      seed=as.integer(seed.id*i))
    
    draws <- Reduce(function(x,y) Map(rbind,x,y), draws.list)
    
    if (any(is.na(draws$counts))) {
        LML <- NA
    } else {
        LML <- get.LML(counts=draws$counts,log.phi=log.phi,
                       post.mode=post.mode,
                       fn.dens.post=get.f.wrap,
                       fn.dens.prop=dmvn.sparse.wrap,
                       prop.params=prop.params)
    }
    ## Section H:  Compute log marginal likelihood
    
    acc.rate <- 1/mean(draws$counts)
    
    dimnames(draws$draws) <- list(iteration=1:NROW(draws$draws),
                                  variable=var.names)
    
    
    draws$LML <- LML
    draws$acc.rate <- acc.rate
    
}
