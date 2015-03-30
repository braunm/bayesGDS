
context("binary")

test_that("small", {

    require(sparseHessianFD)
    require(sparseMVN)
    require(trustOptim)

    seed.id <- 123
    set.seed(seed.id*7)

    rmvn.sparse.wrap <- function(n.draws, params) {
        res <- rmvn.sparse(n.draws, params$mean, params$CH, prec=TRUE)
        return(res)
    }

    dmvn.sparse.wrap <- function(d, params) {
        res <- dmvn.sparse(d, params$mean, params$CH, prec=TRUE)
        return(res)
    }

    data(binary_small)

    D <- binary_small
    N <- length(D[["Y"]])
    k <- NROW(D[["X"]])
    nvars <- as.integer(N*k + k)
    priors <- list(inv.Sigma = rWishart(1,k+5,diag(k))[,,1],
                   inv.Omega = diag(k))
    start <- rnorm(nvars) ## random starting values
    f <- binary.f(start, data=D, priors=priors)

    hs <- drop0(tril(binary.hess(start, data=D, priors=priors)))
    hsNZ <- Matrix.to.Coord(hs)
    FD <- sparseHessianFD.new(start, binary.f, binary.grad,
                              rows=hsNZ$rows, cols=hsNZ$cols,
                              data=D, priors=priors)


?
    f <- FD$fn(start)
    df <- FD$gr(start)
    hess <- FD$hessian(start)

    print("Finding posterior mode")
    opt <- trust.optim(start, fn=FD$fn,
                       gr = FD$gr,
                       hs = FD$hessian,
                       method = "Sparse",
                       control = list(
                           start.trust.radius=5,
                           stop.trust.radius = 1e-7,
                           prec=1e-7,
                           report.precision=1L,
                           maxit=500L,
                           preconditioner=1L,
                           function.scale.factor=-1,
                           report.freq=50L
                       )
                       )

    post.mode <- opt$solution
    hess <- opt$hessian
    var.names <- names(post.mode)

    n.draws <- 300  ## total number of draws needed
    M <- 10000  ## proposal draws
    max.tries <- 100000  ## to keep sample.GDS from running forever
    ds.scale <- .99  ## scaling factor for proposal density

    chol.hess <- Cholesky(-ds.scale*hess)

    prop.params <- list(mean = post.mode,
                        CH = chol.hess
                        )

    log.c1 <- opt$fval
    log.c2 <- dmvn.sparse.wrap(post.mode, prop.params)

    cat("Collecting GDS Proposal Draws\n")
    draws.m <- as(rmvn.sparse.wrap(M,prop.params),"matrix")
    log.post.m <- plyr::aaply(draws.m, 1, FD$fn, .parallel=FALSE)
    log.prop.m <- dmvn.sparse.wrap(draws.m, params=prop.params)
    log.phi <- log.post.m - log.prop.m +log.c2 - log.c1


    invalid.scale <- any(log.phi>0)
    cat("Are any log.phi > 0?  ",invalid.scale,"\n")

    expect_false(invalid.scale)

    if (!invalid.scale) {

        cat("Generating DS draws - accept-reject phase\n")

        draws <- sample.GDS(n.draws = n.draws,
                            log.phi = log.phi,
                            post.mode = post.mode,
                            fn.dens.post = FD$fn,
                            fn.dens.prop = dmvn.sparse.wrap,
                            fn.draw.prop = rmvn.sparse.wrap,
                            prop.params = prop.params,
                            report.freq = 50,
                            thread.id = 1,
                            announce=FALSE)


        expect_false(any(is.na(draws$counts)))


        if (any(is.na(draws$counts))) {
            LML <- NA
        } else {
            LML <- get.LML(counts=draws$counts,
                           log.phi=log.phi,
                           post.mode=post.mode,
                           fn.dens.post= FD$fn,
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


})
