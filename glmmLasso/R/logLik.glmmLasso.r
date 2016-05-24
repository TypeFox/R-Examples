logLik.glmmLasso<-function(y,mu,ranef.logLik=NULL,family,penal=FALSE, K = NULL) 
{
  fam <- family$family

  if(fam=="poisson")
  loglik <- sum(y*log(mu)-mu)
  
  if(fam=="binomial")
  loglik <- sum(log(mu[y==1]))+sum(log((1-mu)[y==0]))
  
  if(fam=="gaussian")
  loglik <- sum(y*mu-0.5*(mu^2))
  
  if(fam=="acat")
  {
  mu_cat <- matrix(mu, byrow = TRUE, ncol = K)
  mu_cat <- cbind(mu_cat,1-rowSums(mu_cat))
  yhelp <- matrix(y, byrow = TRUE, ncol = K)
  yhelp <- cbind(yhelp,1-rowSums(yhelp))
  loglik <- sum(yhelp*log(mu_cat))
  }
  
  if(fam=="cumulative")
  {
    mu_cat <- matrix(mu, byrow = TRUE, ncol = K)  
    mu_help <- mu_cat
      
     
        for (i in 2:K){
          mu_cat[,i] <- mu_help[,i]-mu_help[,i-1]
        }
     
      mu_cat <- cbind(mu_cat,1-mu_help[,K])
      
    yhelp <- matrix(y, byrow = TRUE, ncol = K)
    ynew <- apply(yhelp,1,match,x=1)
    ynew[is.na(ynew)] <- K+1
    ynew <- factor(ynew)
    yhelp <- as.matrix(model.matrix(~0+ynew, contrasts = list(ynew = "contr.treatment")))

      loglik <- sum(yhelp*log(mu_cat))
      loglik
    }

  if(penal)
  loglik <- loglik + ranef.logLik 

return(loglik)
}
