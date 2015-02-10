## deprecated.R -- Part of the bayesGDS package 
## Copyright (C) 2013-2015 Michael Braun
## See LICENSE file for details.



## Place to hold functions that will not be maintained in the future

#' @name bayesGDS-deprecated
#' @aliases vech inv.vech
#' @title Deprecated functions
#' @description These functions were in earlier versions, but will no
#' longer be maintained in this package.  They will likely be moved to
#' another package a some time.
NULL


#' @title vech operator on a square matrix
#' @param M a matrix
#' @return A vector containing the lower triangle of M, ordered column-wise.
#' @rdname bayesGDS-deprecated
#' @export
vech <- function( M )
{
    .Deprecated("matrix")
    
    if ( nrow(M)!=ncol(M))
        stop( "argument M is not a square numeric matrix" )
    return( t( t( M[!upper.tri(M)] ) ) )
}

#' @title inverse vech operator on a vector
#' @param x A vector of conforming length
#' @return A k x k lower triangular matrix
#' @details Returns a lower triangular matrix, with elements
#' determined by x. x must be a vector of length k(k+1)/2, where k is
#' the number of rows (and columns) of the result.
#' @rdname bayesGDS-deprecated
#' @export
inv.vech <- function( x ) {

    .Deprecated("as.vector")
    
    n <- length(x)
    R <- (sqrt(1+8*n)-1)/2

    if (R != as.integer(R)) {
        stop ("in function inv.vech:  vector will not fit in square matrix\n")
    }
    
    res <- matrix(0,R, R)
    res[!upper.tri(res)] <- x
    return(res)
}


#' @title Logit transformation
#' @param p A scalar, vector or matrix, where each element is between
#' 0 and 1.
#' @return result = log(p/(1-p))
#' @rdname bayesGDS-deprecated
#' @export
logit <- function(p) {

  if (any(p<=0) || any(p>=1)) {
    stop (" in function logit:  all elements must be in open interval (0,1)")
  }
  
  res <- log(p) - log(1-p)
  return(res)
  
}

#' @title Inverse logit transformation
#' @param p A scalar, vector or matrix
#' @return result = exp(x)/(1+exp(x))
#' @rdname bayesGDS-deprecated
#' @export
inv.logit <- function(x) {

  ## numerically stable inv.logit
  
  w.max <- x>=log(.Machine$double.xmax)
 
  res <- exp(x - log1p(exp(x)))
  res[w.max] <- 1
  return(res)

}

#' @title Log inverse logit transformation
#' @param p A scalar, vector or matrix
#' @return result = log[exp(x)/(1+exp(x))]
#' @rdname bayesGDS-deprecated
#' @export
log_inv.logit <- function(x) {

  w.max <- x>=log(.Machine$double.xmax)
  w.min <- x<=log(.Machine$double.xmin)
  ww <- !(w.min | w.max)
  
  res <- x
  res[ww] <- x[ww] - log1p(exp(x[ww]))
  res[w.max] <- 0
  return(res)
  
}

