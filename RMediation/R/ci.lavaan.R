ci.lavaan <- function(Obj, quant, alpha=.05, type="MC", plot=TRUE, plotCI=TRUE, n.mc = 1e+06, H0 = FALSE, mu0 = NULL, Sigma0 = NULL){
  fm1 <- Obj
  if(!fm1@Options$do.fit) fm1 <- update(Obj, do.fit=TRUE) # If it's not fitted, refit.
  pEstM1 <- coef(fm1) #parameter estimate
  #quant <- formula(paste("~",FUN))
  name1 <- all.vars(quant)
  coefResM1 <- pEstM1[name1]
  ##Cov 
  covM1 <- (vcov(fm1)) #covariance of the coef estimates
  covResM1 <- covM1[name1,name1] 
  ci(coefResM1, covResM1, quant, alpha, type, plot, plotCI,n.mc,H0, mu0 , Sigma0)
  #AcceptReg <- ci(coefResM0,covResM0,quant,alpha,type,plot,plotCI, add=TRUE)[[1]]
  #return(list(LRT,LRT_Scaled))
}