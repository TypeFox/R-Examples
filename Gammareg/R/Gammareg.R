Gammareg<-function(formula1,formula2,meanlink="log")
{


if (meanlink=="log") {est <- gammahetero1(formula1,formula2)}
else {    
   if(meanlink=="ide") {est <- gammahetero2(formula1,formula2)} 
   else {stop("No link avalaible")}
   }


  est$coefficients <- matrix(c(est$beta,est$gamma))
  names <- c(colnames(est$X),colnames(est$Z))
  par <-rep(c("beta.","gamma."),c(ncol(est$X),ncol(est$Z)))
  rownames(est$coefficients) <-paste(par,names,sep="")
  
  est$desvB <-  est$CovarianceMatrixbeta
  est$desvG <- est$CovarianceMatrixgamma
  est$interv<- rbind(est$ICB,est$ICG)
  est$AIC<-est$AIC  
  est$iteration <-est$iteration  
  est$convergence <-est$convergence

  est$call <- match.call()
  
  class(est) <- "Gammareg"
  
  est 


}
