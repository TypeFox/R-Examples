LD <- function(x, lim, lim2, measure = c("r2", "r", "D"), trim = TRUE) {

  if(!is(x,"bed.matrix")) stop('x must be a bed matrix')
  if(!is.vector(lim))  stop('lim must be a vector of length 2')
  if(length(lim) == 1) lim = c(lim, lim)
  if(length(lim) != 2) stop('lim must be a vector of length 2')
  if(any(lim < 1) | any(lim > ncol(x)) | lim[1] > lim[2]) stop("inappropriate range in lim")

  if(!missing(lim2)) {
    if(!is.vector(lim2))  stop('lim2 must be a vector of length 2')
    if(length(lim2) == 1) lim2 = c(lim2, lim2)
    if(length(lim2) != 2) stop('lim2 must be a vector of length 2')
    if(any(lim2 < 1) | any(lim2 > ncol(x)) | lim2[1] > lim2[2]) stop("inappropriate range in lim2")
  }

  measure <- match.arg(measure) 

  if(!x@standardize_mu_sigma & !x@standardize_p) {
     if(!is.null(x@mu) & !is.null(x@sigma))
       x@standardize_mu_sigma <- TRUE
     else
       stop("Can't center/scale x for LD computation (use set.stats(x))\n")
  }

  if(measure == "D") {
    x@standardize_mu_sigma <- TRUE
    x@sigma <- rep(sqrt(2), ncol(x))
  }  

  if(missing(lim2) || all(lim2 == lim)) {  
    if(x@standardize_mu_sigma) 
      a <- .Call('gg_LD', PACKAGE = 'gaston', x@bed, x@mu, x@sigma, lim[1]-1, lim[2]-1)
    else if(x@standardize_p) {
      warning("Moment estimates of LD using a p-standardized matrix can be outside of the range [-1,1]")
      a <- .Call('gg_LD_p', PACKAGE = 'gaston', x@bed, x@p, lim[1]-1, lim[2]-1)
    }
    rownames(a) <- colnames(a) <- x@snps$id[seq(lim[1], lim[2])]
  } else { 
    if(x@standardize_mu_sigma) 
      a <- .Call('gg_LD_chunk', PACKAGE = 'gaston', x@bed, x@mu, x@sigma, lim[1]-1, lim[2]-1, lim2[1]-1, lim2[2]-1)
    else if(x@standardize_p)
      a <- .Call('gg_LD_chunk_p', PACKAGE = 'gaston', x@bed, x@p, lim[1]-1, lim[2]-1, lim2[1]-1, lim2[2]-1)
    rownames(a) <- x@snps$id[seq(lim[1], lim[2])]
    colnames(a) <- x@snps$id[seq(lim2[1], lim2[2])]
  }

  if(measure != "D" & trim) {
    a[a < -1] <- -1; 
    a[a >  1] <-  1;
  }

  if(measure == "r2") a <- a**2

  a
}

