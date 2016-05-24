"run.old.cor2.examp" <-
function(n=100,seed) {
  if (!missing(seed)){ set.seed(seed) }
  if(!requireNamespace('tcltk', quietly = TRUE)){stop('The tcltk package is needed')}

  x <- scale(matrix(rnorm(2*n,0,1), ncol=2))
  x <- x %*% solve( chol( cor(x) ) )
  xr <- range(x)

  r.old <- 0
  r2.old <- 0

  cor.refresh <- function(...) {
    r <- slider(no=1)
    r2 <- slider(no=2)

    if (r!=r.old){
      slider(set.no.value=c(2,r^2))
      r.old <<- r
      r2.old <<- r^2
    } else {
      slider(set.no.value=c(1, ifelse(r<0, -sqrt(r2), sqrt(r2))))
      r.old <<- ifelse(r<0, -sqrt(r2), sqrt(r2))
      r2.old <<-r2
      r <- r.old
    }

    if ( r == 1 ) {
      cmat <- matrix( c(1,0,1,0),2 )
    } else if (r == -1) {
      cmat <- matrix( c(1,0,-1,0),2 )
    } else {
      cmat <- chol( matrix( c(1,r,r,1),2) )
    }
    new.x <- x %*% cmat

    plot(new.x, xlab='x',ylab='y', xlim=xr, ylim=xr)
    title(paste("r = ",round(cor(new.x[,1],new.x[,2]),3),
                "\nr^2 = ",round(r^2,3)))
  }

  slider( cor.refresh, c('Correlation','r^2'), c(-1,0), c(1,1), c(0.01,0.01),
         c(0,0),
         title="Correlation Demo")
  cor.refresh()
}

