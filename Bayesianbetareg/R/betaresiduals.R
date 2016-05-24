betaresiduals <-
function(Y,X,model){
  
  Y <- as.matrix(Y)
  residuales <- model$residuals
  variance <- model$variance
  phi <- model$precision
  yestimado <-  model$fitted.values
  
  Y.star <- log(Y/(1-Y))
  mu.star <- digamma(yestimado*phi)-digamma((1-yestimado)*phi)
  var.star <- trigamma(yestimado*phi)+trigamma((1-yestimado)*phi)
  
  
  W <- diag(as.numeric(phi*var.star*((yestimado*(1-yestimado)^2))),nrow=nrow(Y))
  H <- (W^0.5)%*%X%*%solve(t(X)%*%W%*%X)%*%t(X)%*%(W^0.5)
  
  #Absolute residuals  
  rabs<-abs(residuales)
  
  #Standardized Weighted Residual 1
  rp<-residuales/sqrt(variance)
  
  #Standardized Weighted Residual 1
  swr1 <- (Y.star - mu.star)/sqrt(var.star)
  
  #Standardized Weighted Residual 2
  swr2 <- (Y.star - mu.star)/sqrt(var.star*(1-diag(H))) 
  
  #Verosimilitud function in t
  lt <- function(Y,phi,mu){
    (dbeta(Y, phi*mu,phi*(1-mu)))
  }
  
  #Residuals deviance
  rtd = sign(residuales)*sqrt(2*abs((lt(Y,phi,Y)-lt(Y,phi,yestimado))))
  
  #Cook's distance
  Cook <- diag(H)*(rp)^2/(ncol(X)*(1-diag(H))^2)   
  
  betaresiduals<- list()
  betaresiduals$abs <- rabs
  betaresiduals$swr0 <-rp
  betaresiduals$swr1 <- swr1
  betaresiduals$swr2 <- swr2
  betaresiduals$deviance <- rtd
  betaresiduals$cook <- Cook
  betaresiduals$H <- H
  
  betaresiduals
  
}
