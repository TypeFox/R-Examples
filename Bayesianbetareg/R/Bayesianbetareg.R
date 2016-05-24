Bayesianbetareg <-
function (Y,X,Z,nsim,bpri,Bpri,gpri,Gpri,burn,jump,graph1,graph2){
  
  est <- BayesianbetaregEst (Y,X,Z,nsim,bpri,Bpri,gpri,Gpri,burn,jump,graph1,graph2)
  
  est$coefficients <- matrix(c(est$Bestimado,est$Gammaest))
  names <- c(colnames(est$X),colnames(est$Z))
  par <-rep(c("beta.","gamma."),c(ncol(est$X),ncol(est$Z)))
  rownames(est$coefficients) <-paste(par,names,sep="")
  
  est$desv <-  matrix(c(est$DesvBeta, est$DesvGamma))
  est$interv<- rbind(est$B,est$G)
  est$fitted.values <- est$yestimado
  est$residuals <- est$residuales
  est$precision <- est$phi
  est$variance <- est$variance
  est$beta.mcmc<-est$beta.mcmc
  est$gamma.mcmc<-est$gamma.mcmc
  est$beta.mcmc.short<-est$beta.mcmc.auto
  est$gamma.mcmc.short<-est$gamma.mcmc.auto
  est$Y<-est$Y
  est$X<-est$X  

  est$call <- match.call()
  
  class(est) <- "Bayesianbetareg"
  
  est 
  
}
