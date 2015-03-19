

##```{r setup1, echo = FALSE}
##knitr::opts_chunk$set(collapse = FALSE, comment = "#", message=FALSE)
##```



##```{r setup2, echo=FALSE}

##----setup2
require(Matrix, quietly=TRUE)
NN <- 6
kk <- 2
pp <- 2
nv1 <- NN*kk+pp
nels1 <- nv1^2
nnz1 <- NN*kk^2 + pp^2 + 2*NN*pp*kk
nnz1LT <- NN*kk*(kk+1)/2 + pp*(pp+1)/2 + pp*NN*kk
Q <- 1000
nv2 <- Q*kk+pp
nels2 <- nv2^2
nnz2 <- Q*kk^2 + pp^2 + 2*Q*kk*pp
nnz2LT <- Q*kk*(kk+1)/2 + pp*(pp+1)/2 + Q*kk*pp


##----pattern1
MM <- as(kronecker(diag(NN),matrix(1,kk,kk)),"lMatrix")
MM <- rBind(MM, Matrix(TRUE,pp,NN*kk))
MM <- cBind(MM, Matrix(TRUE, kk*NN+pp, pp))
print(as(MM,"lgCMatrix"))



##----pattern2
MM <- as(kronecker(matrix(1,kk,kk), diag(NN)),"lMatrix")
MM <- rBind(MM, Matrix(TRUE,pp,NN*kk))
MM <- cBind(MM, Matrix(TRUE, kk*NN+pp, pp))
print(as(MM,"lgCMatrix"))




##----data
##data(binary)
data(binary_test)
binary <- binary_test
str(binary)
N <- length(binary[["Y"]])
k <- NROW(binary[["X"]])
p <- k
nvars <- as.integer(N*k + p)
priors <- list(inv.Sigma = rWishart(1,k+5,diag(k))[,,1],
               inv.Omega = diag(k))



##----startingVals
start <- rnorm(nvars) ## random starting values
f <- binary.f(start, data=binary, priors=priors)
f
df <- binary.grad(start, data=binary, priors=priors)
str(df)
d2f <- binary.hess(start, data=binary, priors=priors)
print(d2f[1:6,1:6], digits=3)


# ---- hessStruct
require(sparseHessianFD, quietly=TRUE)
hs <- Matrix(0, nvars, nvars)
for (i in 1:(N + 1)) {
    ## range of row / col indices of block diagonal
    rng <- ((i-1)*k+1):(k*i)
    hs[rng, rng] <- tril(Matrix(1,k,k)) ## lower triangle
}
hs[N*k + 1:p, 1:(N*k)] <- 1 ## bottom margin
hsNZ <- Matrix.to.Coord(hs)
str(hsNZ)


##----hessStruct2
## require(sparseHessianFD)
## hs <- drop0(tril(binary.hess(start, data=binary, priors=priors)))
## hsNZ <- Matrix.to.Coord(hs)
## str(hsNZ)

##----sparseHessianFD
FD <- sparseHessianFD.new(start, binary.f, binary.grad,
                          rows=hsNZ[["rows"]], cols=hsNZ[["cols"]],
                          data=binary, priors=priors)



##----usingFD
f <- FD$fn(start)
df <- FD$gr(start)
hess <- FD$hessian(start)



##----hessUpperLeft
print(hess[1:6,1:6], digits=3)
all.equal(hess, d2f)


##----trustOptim
require(trustOptim, quietly=TRUE)
opt <- trust.optim(start, fn=FD$fn, gr = FD$gr, hs = FD$hessian,
                   method = "Sparse",
                   control = list(
                       start.trust.radius=5, stop.trust.radius = 1e-7,
                       prec=1e-7, report.precision=1,
                       maxit=500, preconditioner=1,
                       function.scale.factor=-1
                       )
                   )

theta.star<- opt[["solution"]]
hess <- opt[["hessian"]]
var.names <- names(theta.star)



##----defPropFuncs
require(sparseMVN, quietly=TRUE)
rmvn.sparse.wrap <- function(n.draws, params) {
## sample MVN with sparse precision matrix
    res <- rmvn.sparse(n.draws, params[["mean"]], params[["CH"]], prec=TRUE)
    return(res)
}

dmvn.sparse.wrap <- function(d, params) {
## MVN density with sparse precision
    res <- dmvn.sparse(d, params[["mean"]], params[["CH"]], prec=TRUE)
    return(res)
}



##----propParams
scale <- .96
chol.hess <- Cholesky(-scale*hess)
prop.params <- list(mean = theta.star, CH = chol.hess)

##----parallelSetup
library(doParallel, quietly=TRUE)
run.par <- TRUE
if(run.par) registerDoParallel(cores=10) else registerDoParallel(cores=1)
seed.id <- 123
set.seed(seed.id)



##----proposals
M <- 10000  ## proposal draws
log.c1 <- FD$fn(theta.star)
log.c2 <- dmvn.sparse.wrap(theta.star, prop.params)
draws.m <- as(rmvn.sparse.wrap(M,prop.params),"matrix")
log.post.m <- plyr::aaply(draws.m, 1, FD$fn, .parallel=run.par)
log.prop.m <- dmvn.sparse.wrap(draws.m, params=prop.params)
log.phi <- log.post.m - log.prop.m + log.c2 - log.c1
valid.scale <- all(log.phi <= 0)
cat("Are all log.phi <= 0?  ",valid.scale,"\n")




##----sampleGDS
## if valid.scale is FALSE, need to change
## the proposal density

n.draws <- 20  ## total number of draws needed
max.tries <- 100000  ## to keep sample.GDS from running forever
if (valid.scale) {
   if (run.par) {
   ## running in parallel, 1 sample per core
      batch.size <- 2
      n.batch <- floor(n.draws / batch.size)
      draws.list <- foreach(i=1:n.batch, .inorder=FALSE) %dopar% sample.GDS(
                                               n.draws = n.draws,
                                               log.phi = log.phi,
                                               post.mode = theta.star,
                                               fn.dens.post = FD$fn,
                                               fn.dens.prop = dmvn.sparse.wrap,
                                               fn.draw.prop = rmvn.sparse.wrap,
                                               prop.params = prop.params,
                                               report.freq = 50,
                                               thread.id = i,
                                               announce=TRUE,
                                               seed=as.integer(seed.id*i))

        ## combine results from each batch
        draws <- Reduce(function(x,y) Map(rbind,x,y), draws.list)
    } else {
        ## single core
        draws <- sample.GDS(n.draws = n.draws, log.phi = log.phi,
                            theta.star = theta.star,   fn.dens.post = FD$fn,
                            fn.dens.prop = dmvn.sparse.wrap,
                            fn.draw.prop = rmvn.sparse.wrap,
                            prop.params = prop.params,
                            report.freq = 50, thread.id = 1, announce=TRUE)
    }
}


##----strDraws
str(draws)


##----summary
quants <-  plyr::aaply(draws[["draws"]][,(N*k+1):NCOL(draws[["draws"]])], 2,
       quantile, probs=c(.025, .5, .975))
quants




##----LML
if (valid.scale) {
    if (any(is.na(draws[["counts"]]))) {
        LML <- NA
    } else {
        LML <- get.LML(counts=draws$counts,
                       log.phi=log.phi,
                       post.mode=theta.star,
                       fn.dens.post= FD$fn,
                       fn.dens.prop=dmvn.sparse.wrap,
                       prop.params=prop.params)
    }
    acc.rate <- 1/mean(draws$counts)
    dimnames(draws$draws) <- list(iteration=1:NROW(draws[["draws"]]),
                                  variable=var.names)
    draws[["LML"]] <- LML
    draws[["acc.rate"]] <- acc.rate
}
