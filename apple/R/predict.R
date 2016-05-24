predict.apple <- function(object,X,which=1:length(object$lambda), type=c("link","response","class"), ...)
  {
    #if (type=="coefficients") return(object$beta[,which])
    eta <- t(object$a0[which] + t(X%*%object$beta[,which]))
    if (type=="link") return(eta)
    if(object$family=='binomial'){
    	pihat <- exp(eta)/(1+exp(eta))
        }
    if(object$family=='poisson'){
    	pihat <- exp(eta)
        }
    if (type=="response") return(pihat)
    if(type=='class' & object$family=='binomial') return(eta>0)
    if(type=='class' & object$family=='poisson') stop('choose response for Poisson, or choose binomial for response')
  }
