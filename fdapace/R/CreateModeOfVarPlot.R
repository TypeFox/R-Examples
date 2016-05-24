#' Functional Principal Component Analysis mode of variation plot
#' 
#' Create the k-th mode of variation plot around the mean. The red-line is
#' the functional mean, the grey shaded areas show the range of variations
#' around the mean: \eqn{ \pm Q \sqrt{\lambda_k} \phi_k}{+/- Q sqrt{lambda_k} phi_k}
#' for the dark grey area Q = 1, and for the light grey are Q = 2.
#'
#' @param fpcaObj An FPCA class object returned by FPCA(). 
#' @param k The k-th mode of variation to plot (default k = 1) 
#' @param ... Additional arguments for the 'plot' function.
#'
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- Wiener(n, pts)
#' sampWiener <- Sparsify(sampWiener, pts, 10)
#' res <- FPCA(sampWiener$yList, sampWiener$tList, 
#'             list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
#' CreateModeOfVarPlot(res)
#' @export

CreateModeOfVarPlot <-function(fpcaObj,  k = 1, ...){ 
  
  args1 <- list( main="Default Title", xlab='s', ylab='')  
  inargs <- list(...)
  args1[names(inargs)] <- inargs
  
  if(k> length(fpcaObj$lambda) ){
    stop("You are asking to plot a mode of variation that is incomputable.")
  }  
  
  obsGrid = fpcaObj$obsGrid      
  s = fpcaObj$workGrid
  mu = fpcaObj$mu
  
  sigma = sqrt(fpcaObj$lambda[k])
  sigma1 = sqrt(fpcaObj$lambda[1])
  phi = fpcaObj$phi[,k]
  phi1 = fpcaObj$phi[,1]
  
  do.call(plot, c(list(type='n'), list(x=s), list(y=s), 
                  list(ylim=range(c( 3* sigma1 * phi1 + mu , -3* sigma1 * phi1 + mu ))), args1))
  grid()    
  polygon(x=c(s, rev(s)), y = c( -2* sigma * phi + mu, 
                                 rev(2* sigma * phi + mu)), col= 'lightgrey',border=NA)
  polygon(x=c(s, rev(s)), y = c( -1* sigma * phi + mu, 
                                 rev(1* sigma * phi + mu)), col= 'darkgrey',border=NA)  
  lines(x=s, y=mu , col='red')
}
