manipulate.cor.examp <- function(n=100, 
                                 seed) {
  if(!requireNamespace('manipulate', quietly = TRUE)) stop("This function depends on the manipulate package within Rstudio")
  
  if(!missing(seed)) set.seed(seed)
  
  x <- scale(matrix(rnorm(2*n,0,1), ncol=2))
  x <- x %*% solve( chol( cor(x)))
  xr <- range(x,-x)
  
  replot <- function(r) {
    if( r >= 1 ) {
      cmat <- matrix( c(1,0,1,0), 2 )
    } else if( r <= -1 ) {
      cmat <- matrix( c(1,0,-1,0), 2 )
    } else {
      cmat <- chol( matrix(c(1,r,r,1),2) )
    }
    new.x <- x %*% cmat
    plot(new.x, xlab='x',ylab='y', 
         xlim=xr,ylim=xr)
    title(paste("r =", round(r,3)))
  }
  
  r <- NA # so that following function does not complain about global var
  manipulate::manipulate(replot(r),r=manipulate::slider(-1,1,0,step=0.005))
  
}