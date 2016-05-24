closedp.Mtb <- function(X, dfreq = FALSE, method = "BFGS", ...)
{
  
  ############################################
  # Validation des arguments fournis en entrée
  valid.one(dfreq,"logical")
  Xvalid<-valid.X(X,dfreq)
    X <- Xvalid$X
    t <- Xvalid$t
  if (t < 3) stop("model Mtb cannot be fitted with 2 capture occasions")  
  ############################################

  Y <- histfreq.t(X,dfreq=dfreq)
  if (any(Y>170)) stop("model Mtb cannot be fitted with capture history frequencies larger than 170")  
  
  n <- sum(na.rm=TRUE,Y)
  histpos <- histpos.t(t)
  nbcap <- rowSums(histpos)
        
  pnames <- paste("p",1:t,sep="")
  cnames <- paste("c",2:t,sep="")


  ### Modele Mtb
  
  # matrices d'indicatices de non capture (wi=0) avant la première capture, i=1,...,t-1
  I1 <- matrix(c(0,0,1),ncol=1)  
  for(i in (3:t))
    I1 <- cbind(c(rep(0, 2^(i-1)), rep(1,((2^(i-1))-1))), 
                rbind(matrix(rep(0,(2^(i-1))*(i-2)),nrow=2^(i-1)),I1))
  # matrice d'indicatrices de la première capture
  I2 <- matrix(rep(0,(2^t-1)*t),ncol=t)  
  poscapt1<-rep(1,2^(t-1))
  for(i in 2:t) poscapt1 <- c(poscapt1,rep(i,2^(t-i)))
  for(i in 1:(2^t-1)) I2[i,poscapt1[i]] <- 1        
  # matrice d'indicatrices de capture (wi=1) après la première capture
  I3 <- histpos.t(t)[,-1]  
  for(i in (2^(t-1)+1):(2^t-1)) 
    I3[i,1:(poscapt1[i]-1)]<-rep(0,length(I3[i,1:(poscapt1[i]-1)]))
  # matrice d'indicatrices de non-capture (wi=0) après la première capture
  I4 <- 1-histpos.t(t)[,-1]  
  for(i in (2^(t-1)+1):(2^t-1)) 
    I4[i,1:(poscapt1[i]-1)]<-rep(0,length(I4[i,1:(poscapt1[i]-1)]))
  
  # valeurs initiales des paramètres (N->n, tous les p -> 0.5 et ci=pi)
  init_loglinparam.Mtb <- as.vector(c(log(n),rep(0,t+1)))
  
  # Fonction de log-vraisemblance
  ll_loglinparam.Mtb <- function(par,Y,I1,I2,I3,I4)
  {
          logmu <- par[1] + I1%*%log(1/(1+exp(par[2:t]))) + 
            I2%*%log(exp(par[2:(t+1)])/(1+exp(par[2:(t+1)]))) + 
            I3%*%log(exp(par[3:(t+1)]+par[t+2])/(1+exp(par[3:(t+1)]+par[t+2]))) + 
            I4%*%log(1/(1+exp(par[3:(t+1)]+par[t+2])))
          # Sortie : log-vraisemblance
          t(Y)%*%logmu - sum(na.rm=TRUE,exp(logmu)) - sum(na.rm=TRUE,log(factorial(Y))) 
  }
  
  # Optimisation de la fonction de log-vraisemblance en fonction des paramètres loglinéaires
  optim.out <- optim.call(par = init_loglinparam.Mtb, fn = ll_loglinparam.Mtb, method=method,
                          Y = Y, I1 = I1, I2 = I2, I3 = I3, I4 = I4, 
                          control=list(fnscale = -1, maxit = 1000), ...)
  
  # Si la commande à généré une erreur, on arrête l'exécution de la fonction
  if (!is.null(optim.out$error))
    stop("error when fitting the model: ", optim.out$error)  
  
  # Statistiques d'ajustement du modèle
  optimo <- optim.out$optimo
  logmu <- optimo$par[1] + I1%*%log(1/(1+exp(optimo$par[2:t]))) + 
    I2%*%log(exp(optimo$par[2:(t+1)])/(1+exp(optimo$par[2:(t+1)]))) + 
    I3%*%log(exp(optimo$par[3:(t+1)]+optimo$par[t+2])/(1+exp(optimo$par[3:(t+1)]+optimo$par[t+2]))) + 
    I4%*%log(1/(1+exp(optimo$par[3:(t+1)]+optimo$par[t+2])))
  se <- sqrt(exp(2*optimo$par[1])*solve(-optimo$hessian)[1,1]-exp(optimo$par[1]))
  dev <- 2*sum(na.rm=TRUE,Y*(pmax(log(Y),0)-logmu))
  df <- 2^t-1-(t+2)
  AIC <- -2*optimo$value + 2*(t+2)
  BIC <- -2*optimo$value + log(n)*(t+2)
  
  # Préparation des sorties
  infoFit <- getInfo(err = optim.out$error, warn = optim.out$warnings)
  tableau <- matrix(c(exp(optimo$par[1]),se,dev,df,AIC,BIC,infoFit),nrow=1)
  
  # Autres paramètres
  parMtb <- matrix(c(exp(optimo$par[1]), exp(optimo$par[2:(t+1)])/(1+exp(optimo$par[2:(t+1)])),
                     exp(optimo$par[3:(t+1)]+optimo$par[t+2])/
                       (1+exp(optimo$par[3:(t+1)]+optimo$par[t+2]))),nrow=1)
  colnames(parMtb) <- c("N",pnames,cnames)


  # Préparation des sorties
  dimnames(tableau) <- list("Mtb",c("abundance","stderr","deviance","df","AIC","BIC","infoFit"))
  ans <- list(n=n, t=t, results=tableau, optim = optimo, optim.warn = optim.out$warnings, 
              parMtb=parMtb)
  class(ans) <- "closedp.Mtb"
  ans

}


print.closedp.Mtb <- function(x, ...) {
  cat("\nNumber of captured units:",x$n,"\n\n")
  cat("Abundance estimation and model fit:\n")
  tabprint(tab = x$results, digits = c(1,1,3,0,3,3,NA), warn = x$optim.warn, ...)
  cat("\nNote: Model Mtb is not loglinear.\nThe abundance estimation for this model can be unstable.\n")
  cat("\n")
  invisible(x)
}
