
get.cutoffs <- function(log.phi, n.draws) {

    ## INPUTS:
    ##   log.phi - log of \Phi (should all be negative)
    ##   n.draws - number of draws required from the posterior
    
    
    ## OUTPUTS:
    ##   log.u - threshold draws, using the log transformation in Braun and Damien, for use in GDS A-R step
    
    if (any(log.phi>0) ) stop("all values of log phi must be non-positive")
    
    M <- length(log.phi)
    v <- c(-log.phi[order(log.phi,decreasing=TRUE)],Inf)
    v1 <- min(v)
    
    p.int <- as.numeric(exp(log(seq(1,M)) -v1 + log(exp(v1-v[1:M]) - exp(v1-v[2:(M+1)]))))
    
    j <- sample(1:M,n.draws,replace=TRUE,prob=p.int)
    
    u <- runif(n.draws)
    d <- as.numeric(v[j] - log(1-u + u*exp(-v[j+1]))) ## draw from right-censored exponential
    
    return(d)
}

.sample.GDS.unc <- function(n.draws,
                           log.phi,
                           post.mode,
                           fn.dens.post,
                           fn.dens.prop,
                           fn.draw.prop,
                           prop.params,
                           ...,
                           max.tries=1000000,        
                           report.freq=1,
                           announce=FALSE,
                           thread.id=1,
                           seed=.Random.seed
                           ) {
   

    if (any(!is.finite(seed))) {
        stop("Error in sample.GDS:  all values in seed must be finite")
    }
    set.seed(seed)
    
    v <- get.cutoffs(log.phi, n.draws)
    nvars <- length(post.mode)
    log.c1 <- fn.dens.post(post.mode, ...)
    log.c2 <- fn.dens.prop(post.mode, prop.params)
    
    v.keep <- vector("numeric",length=n.draws)
    counts <- vector("integer",length=n.draws)
    gt.1 <- vector("integer",length=n.draws)
    draws.keep <- matrix(nrow=n.draws, ncol=nvars)
    log.post.dens <- vector("numeric",length=n.draws)
    log.prop.dens <- vector("numeric",length=n.draws)
    crit.keep <- vector("numeric",length=n.draws)
    
    remaining.draws <- n.draws
    idx.1 <- 1
    count.idx <- 1

    while((remaining.draws>0) & (count.idx <= max.tries)) {

        if ((count.idx %% report.freq) == 0) {
            cat("thread ",thread.id,"  count ",count.idx,"  remaining draws = ",remaining.draws,"\n")
        }
        
        x <- fn.draw.prop(remaining.draws, prop.params)
        dens.1 <- as.vector(apply(x,1,fn.dens.post,...))     
        dens.2 <- as.vector(fn.dens.prop(x,prop.params))
        crit <- dens.1-dens.2-log.c1+log.c2
        ww <- which(-v < crit)
        n.keep <- length(ww)
        if (any(crit>0)) {
            cat("bayesGDS warning:  accepted draw(s) with log.phi>0\n")
        }
        
        if (n.keep > 0) {
            idx.range <- idx.1:(idx.1 + n.keep - 1)
            v.keep[idx.range] <- v[ww]
            counts[idx.range] <- count.idx
            draws.keep[idx.range,] <- x[ww,]
            log.post.dens[idx.range] <- dens.1[ww]
            log.prop.dens[idx.range] <- dens.2[ww]
            gt.1[idx.range] <- crit[ww]>0
            crit.keep[idx.range] <- crit[ww]    
            idx.1 <- idx.1+n.keep
            remaining.draws <- remaining.draws-n.keep
            v <- v[-ww] ## drop v that were successful
            if (announce) {
                cat("Sample accepted in thread ",thread.id,". count =  ",count.idx,"\n")
            }
        }
        count.idx <- count.idx+1
    }
    
    res <- list(draws=draws.keep,
                counts=counts,
                gt.1=gt.1,
                log.post.dens=log.post.dens,
                log.prop.dens=log.prop.dens,
                log.thresholds=v.keep,
                log.phi=crit.keep)
    
    return(res)
}

sample.GDS <- cmpfun(.sample.GDS.unc)
