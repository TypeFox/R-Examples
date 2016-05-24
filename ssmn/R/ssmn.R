ssmn <- function(y, X, family="sn", method="EM", error =  1e-6, maxit=1000, show.envelope = FALSE)
{
  #Running the algorithm
  if (method=="EM"){
  out <- ssmn.est(y, X, family, method, error, maxit)
  cat('\n')
  cat('------------------------------------------------------------\n')
  cat('Skew Scale Mixtures of Normal Distributions\n')
  cat('------------------------------------------------------------\n')
  cat('\n')
  cat('Observations =',length(y))
  cat('\n')
  cat('\n')
  cat('Family =',family)
  cat('\n')
  cat('-------------------------\n')
  cat('Estimates via EM algoritm\n')
  cat('-------------------------\n')
  cat('\n')
  print(round(out$ttable,5))
  cat('\n')

  cat('------------------------\n')
  cat('Model selection criteria\n')
  cat('------------------------\n')
  cat('\n')
  critFin <- c(out$loglik, out$aic, out$bic)
  critFin <- round(t(as.matrix(critFin)),digits=3)
  dimnames(critFin) <- list(c("Value"),c("Loglik", "AIC", "BIC"))
  print(critFin)
  cat('\n')

  cat('-------\n')
  cat('Details\n')
  cat('-------\n')
  cat('\n')
  cat("Convergence reached? =",(out$iter < maxit))
  cat('\n')
  cat('EM iterations =',out$iter,"/",maxit)
  cat('\n')
  cat('Criteria =',round(out$criterio,9))
  cat('\n')
  cat("Processing time =",out$time,units(out$time))
  cat('\n\n')
  }

  if (method=="MLE"){
    out <- ssmn.est(y, X, family, method, error, maxit)
    cat('\n')
    cat('------------------------------------------------------------\n')
    cat('Skew Scale Mixtures of Normal Distributions\n')
    cat('------------------------------------------------------------\n')
    cat('\n')
    cat('Observations =',length(y))
    cat('\n')
    cat('\n')
    cat('Family =',family)
    cat('\n')
    cat('-------------------------------------\n')
    cat('MLE estimates via quasi-Newton method\n')
    cat('-------------------------------------\n')
    cat('\n')
    print(round(out$ttable,5))
    cat('\n')

    cat('------------------------\n')
    cat('Model selection criteria\n')
    cat('------------------------\n')
    cat('\n')
    critFin <- c(out$loglik, out$aic, out$bic)
    critFin <- round(t(as.matrix(critFin)),digits=3)
    dimnames(critFin) <- list(c("Value"),c("Loglik", "AIC", "BIC"))
    print(critFin)
    cat('\n')

    cat('-------\n')
    cat('Details\n')
    cat('-------\n')
    cat('\n')
    cat("Convergence reached? =",(out$convergence))
    cat('\n')
    cat('MLE iterations =',out$iter)
    cat('\n')
    cat('Criteria =',round(out$criterio,9))
    cat('\n')
    cat("Processing time =",out$time,units(out$time))
    cat('\n\n')
  }

  return(out)
}




#########################################################
#Teste

##setwd("~/Dropbox/Library_SSMN/ssmn/R")
#source("functions.R")
#source("ssmn.R")
#source("ssmn.est.R")


#family = "scn" #  "sn"   "stn"  "scn"  "sep" "ssl"
#n0=100
#y = rssmn(n=n0, location=0, scale=1, shape=3, nu=0.6, gama=0.5, dp=NULL, family)
#y = y$y
#X = cbind(rep(1,n0))
#objetos = ssmn(y, X, family="sep", method="EM", error =  1e-6, maxit=5000)
#objetos = ssmn(y, X, family="sep", method="MLE", error =  1e-6, maxit=1000)









