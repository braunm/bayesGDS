## utilities to support generalized direct sampling

vech <- function( M )
{
    if ( nrow(M)!=ncol(M))
        stop( "argument M is not a square numeric matrix" )
    return( t( t( M[!upper.tri(M)] ) ) )
}


inv.vech <- function( x ) {

  n <- length(x)
  R <- (sqrt(1+8*n)-1)/2

  if (R != as.integer(R)) {
    stop ("in function inv.vech:  vector will not fit in square matrix\n")
  }

  res <- matrix(0,R, R)
  res[!upper.tri(res)] <- x
  return(res)
}


logit <- function(p) {

  if (any(p<=0) || any(p>=1)) {
    stop (" in function logit:  all elements must be in open interval (0,1)")
  }
  
  res <- log(p) - log(1-p)
  return(res)
  
}

inv.logit <- function(x) {

  ## numerically stable inv.logit
  
  w.max <- x>=log(.Machine$double.xmax)
 
  res <- exp(x - log1p(exp(x)))
  res[w.max] <- 1
  return(res)

}

log_inv.logit <- function(x) {

  w.max <- x>=log(.Machine$double.xmax)
  w.min <- x<=log(.Machine$double.xmin)
  ww <- !(w.min | w.max)
  
  res <- x
  res[ww] <- x[ww] - log1p(exp(x[ww]))
  res[w.max] <- 0
  return(res)
  
}

