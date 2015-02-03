
get.LML <- function(counts, log.phi, post.mode, fn.dens.post, fn.dens.prop, prop.params,
                    ...) {

  if ( any(!is.finite(counts)) | any(counts<=0) ) {
    stop ("Error in get.LML:  all counts must be finite and positive")
  }

  if ( any(!is.finite(log.phi)) | any(log.phi>0) ) {
    stop ("Error in get.LML:  all values of log.phi must be non-positive")
  }

  log.c1 <- fn.dens.post(post.mode, ...)
  log.c2 <- fn.dens.prop(post.mode, prop.params)
  ar <- 1/mean(counts)
  M <- length(log.phi)
  ord <- order(log.phi,decreasing=TRUE)
  v <- log.phi[ord]
  max.v <- max(v)
  z <- sum((2*seq(1,M)-1) * exp(v-max.v))

  LL <- max.v + log.c2 - log.c1 - log(ar) - 2*log(M) + log(z)
  return(as.numeric(LL))
 

}

