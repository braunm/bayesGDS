library(sparseHessianFD)
library(trustOptim)
library(sparseMVN)
# library(bayesGDS)
library(devtools)
load_all('/Users/rkalescky/iCloud/Source Code/Repositories/bayesGDS')

data(binary_small)
binary <- binary_small
str(binary)
N <- length(binary[["Y"]])
k <- NROW(binary[["X"]])
q <- k
nvars <- as.integer(N*k + q)
priors <- list(inv.Sigma = diag(k), inv.Omega = diag(k))
data.frame(parameter = c("N","k","q"), value = c(N, k, q))

start <- rnorm(nvars) ## random starting values
f <- binary.f(start, data=binary, priors=priors)
df <- binary.grad(start, data=binary, priors=priors)
str(df)

d2f <- binary.hess(start, data=binary, priors=priors)
print(d2f[1:6,1:6], digits=3)

hs <- Matrix(0, nvars, nvars)
for (i in 1:(N + 1)) {
   rng <- ((i-1)*k+1):(k*i)
   hs[rng, rng] <- tril(Matrix(1,k,k))
}
hs[N*k + 1:q, 1:(N*k)] <- 1
hsNZ <- Matrix.to.Coord(hs)
str(hsNZ)

FD <- sparseHessianFD(start, binary.f, binary.grad, rows=hsNZ[["rows"]], cols=hsNZ[["cols"]], data=binary, priors=priors)

f <- FD$fn(start)
df <- FD$gr(start)
hess <- FD$hessian(start)

print(hess[1:6,1:6], digits=3)
all.equal(hess, d2f, tolerance = 1e-7)

opt <- trust.optim(start, fn=FD$fn, gr = FD$gr, hs = FD$hessian,
                   method = "Sparse",
                   control = list(start.trust.radius=5, stop.trust.radius = 1e-7,
                                  prec=1e-7, report.precision=1,
                                  maxit=500, preconditioner=1,
                                  function.scale.factor=-1))

theta.star <- opt[["solution"]]
hess <- opt[["hessian"]]

rmvn.sparse.wrap <- function(n.draws, params) {
    rmvn.sparse(n.draws, params[["mean"]], params[["CH"]], prec=TRUE)
}
dmvn.sparse.wrap <- function(d, params) {
    dmvn.sparse(d, params[["mean"]], params[["CH"]], prec=TRUE)
}

scale <- .96
chol.hess <- Cholesky(-scale*hess)
prop.params <- list(mean = theta.star, CH = chol.hess)

seed.id <- 123
set.seed(seed.id)

M <- 10000
log.c1 <- FD$fn(theta.star)
log.c2 <- dmvn.sparse.wrap(theta.star, prop.params)
draws.m <- as(rmvn.sparse.wrap(M,prop.params),"matrix")
log.post.m <- plyr::aaply(draws.m, 1, FD$fn, .parallel=TRUE)
log.prop.m <- dmvn.sparse.wrap(draws.m, params=prop.params)
log.phi <- log.post.m - log.prop.m + log.c2 - log.c1
valid.scale <- all(log.phi <= 0)
stopifnot(valid.scale)

n.draws <- 1000
max.tries <- 1000000
# if (!run.par) {
#     draws <- sample.GDS(n.draws = n.draws,
#         log.phi = log.phi,
#         theta.star = theta.star,
#         fn.dens.post = FD$fn,
#         fn.dens.prop = dmvn.sparse.wrap,
#         fn.draw.prop = rmvn.sparse.wrap,
#         prop.params = prop.params,
#         report.freq = 1)
# }
# if (run.par) {
#     n.batch <- 10
#     batch.size <- ceiling(n.draws / n.batch)
#     draws.list <- foreach(i=1:n.batch, .inorder=FALSE) %dopar% sample.GDS(
#         n.draws = n.draws,
#         log.phi = log.phi,
#         post.mode = theta.star,
#         fn.dens.post = FD$fn,
#         fn.dens.prop = dmvn.sparse.wrap,
#         fn.draw.prop = rmvn.sparse.wrap,
#         prop.params = prop.params,
#         report.freq = 1,
#         thread.id = i,
#         seed=as.integer(seed.id*i))
#     draws <- Reduce(function(x,y) Map(rbind,x,y), draws.list)
# }

draws <- sample.GDS.rdsm(n.draws = n.draws,
    max.tries = max.tries,
    log.phi = log.phi,
    post.mode = theta.star,
    fn.dens.post = FD$fn,
    fn.dens.prop = dmvn.sparse.wrap,
    fn.draw.prop = rmvn.sparse.wrap,
    prop.params = prop.params,
    report.freq = 1,
    announce = FALSE,
    nodes = 1,
    threadspernode = 1)

str(draws)

quants <-  plyr::aaply(draws[["draws"]][,(N*k+1):nvars], 2, quantile, probs=c(.025, .5, .975), .parallel = FALSE)

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

draws[["LML"]] <- LML
draws[["acc.rate"]] <- 1/mean(draws$counts)
cat("Acceptance rate: ",draws[["acc.rate"]], "\nLog marginal likelihood: ",draws[["LML"]],"\n")

#quit()

