## fast approximation to the full Gamma MLE via Minka (2002) at MS Research
gammaMLE <- function(x,w=NULL,niter=100,tol=0.1,minx=1) { 

  if( is.null(w) ) w <- rep( 1, length(x) )
  meanlogx <- weighted.mean(log(pmax(x,minx)), w)
  meanx <- weighted.mean(pmax(x,minx), w)
  logmeanx <- log(meanx)
  a <- a0 <- (0.5/(logmeanx-meanlogx))  # from Minka 2002
  if(is.nan(a)) stop('NaN starting estimate')
  update.a <- function(a) {
    ooa <- 1/a
    1/(ooa+((meanlogx-logmeanx+log(a)-digamma(a))/(((ooa-trigamma(a))*(a**2)))))
  }
  for(i in 1:niter) { # usually converges in under 5 iterations
    a <- abs(update.a(a0))
    if(abs(a0-a) < tol) break
    else a0 <- a 
  }
  b <- meanx/a
  return(c(shape=a, scale=b))
} 
