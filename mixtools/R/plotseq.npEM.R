# plots SEQuence from npEM object
# function for plotting the scalar parameters sequences along iterations
# if object x comes from npEM, just the x$lambda sequence
# if object x comes from spEM, both x$lambda and x$mu sequences
plotseq <- function(x, ...) UseMethod("plotseq") 



plotseq.npEM <- function(x, ...) {
#  ask <- par(ask=TRUE)
  r <- NCOL(x$data)
  n <- NROW(x$data)
  m <- length(x$lambdahat)
  iter <- NROW(x$lambda)
  xlabel <- paste("iterations")
  nbcol <- 1
  if (!is.null(x$symmetric) && x$symmetric)
		nbcol <- 2
	par(mfcol=c(m, nbcol))

# in all cases, plots the lambda's  
  for (j in 1:m) {
  	estim <- paste(round(x$lambdahat[j],3))
  	tt <- substitute(expression(paste("sequence of ",lambda[j],
  			", estimate ",widehat(lambda[j]),"=", estim, sep="")))
  	ylabel <- substitute(expression(paste(lambda[j],sep="")))
  	plot(x$lambda[,j], type="l", main=eval(tt), xlab=xlabel, 
  		ylab=eval(ylabel), ...)
  	lines(c(0,iter),rep(x$lambdahat[j],2),col=2,lty=2)
  	}

## for symmetric location spEM case plots mu
    if (!is.null(x$symmetric) && x$symmetric) {
  	for (j in 1:m) {
  		estim <- paste(round(x$muhat[j],3))
  	  	tt <- substitute(expression(paste("sequence of ",mu[j],
  	  		", estimate ",widehat(mu[j]),"=",estim,sep="")))
  		ylabel <- substitute(expression(paste(mu[j],sep="")))

  	plot(x$mu[,j], type="l", main=eval(tt),  ylab=eval(ylabel),
  			 xlab=xlabel, ...)
  	lines(c(0,iter),rep(x$muhat[j],2),col=2,lty=2)
   		}
#   	legend("topright", legend=round(x$muhat,3), fill=2:(m+1)) 
    } 
#      structure (list (call=match.call()))
}


