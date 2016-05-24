"run.old.cor.examp" <-
function(n=100,seed) {
  if (!missing(seed)){ set.seed(seed) }
  if(!requireNamespace('tcltk', quietly = TRUE)){stop('The tcltk package is needed')}

  x <- scale(matrix(rnorm(2*n,0,1), ncol=2))
  x <- x %*% solve( chol( cor(x) ) )
  xr <- range(x)

  cor.refresh <- function(...) {
    r <- slider(no=1)
    if ( r == 1 ) {
      cmat <- matrix( c(1,0,1,0),2 )
    } else if (r == -1) {
      cmat <- matrix( c(1,0,-1,0),2 )
    } else {
      cmat <- chol( matrix( c(1,r,r,1),2) )
    }
    new.x <- x %*% cmat

    plot(new.x, xlab='x',ylab='y', xlim=xr, ylim=xr)
    title(paste("r = ",round(cor(new.x[,1],new.x[,2]),3)))
  }

  slider( cor.refresh, 'Correlation', -1, 1, 0.01, 0,
         title="Correlation Demo")
  cor.refresh()
}

