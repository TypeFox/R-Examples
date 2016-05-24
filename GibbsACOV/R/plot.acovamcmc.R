#' @param x A \code{acovamcmc} object
#' @param ... Ignored
#' @method plot acovamcmc

plot.acovamcmc <- function(x, ...){
  #variance plots: density, autocorrelation
  dev.new()
  par(mfrow=c(2,1))
  plot(density(x$sig2a),main="Density of Sigma_a^2")
  plot(density(x$sig2e),main="Density of Sigma_e^2")
  dev.new()
  par(mfrow=c(2,1))
  acf(x$sig2a,main="Series of Sigma_a^2")
  acf(x$sig2e,main="Series of Sigma_e^2")
  #beta plots: auto correlation
  nbetas = nrow(x$Credible_Interval)-2
  dev.new()
  par(mfrow=c(nbetas,1))
  x$beta=as.matrix(x$beta)
  for(i in 1:nbetas){
    acf(x$beta[,i],main=paste("Series of Beta", i))
  }
  
  #beta plots if nbeta = 2
  if(nbetas==2){
  dev.new()
  beta.kde=kde2d(x$beta[,1],x$beta[,2],n=50)
  par(mfrow=c(2,2))
  contour(beta.kde)
  image(beta.kde)
  persp(beta.kde,phi=45,theta=25)
  persp(beta.kde,phi=45,theta=30,shade=0.1,border=NA)
  }
}
