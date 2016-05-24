##' mcmcPriors function
##'
##' A function to define priors for the MCMC.
##' 
##' The package \code{spatsurv} only provides functionality for the built-in Gaussian priors. However, the choice of prior is 
##' extensible by the user by creating functions similar to the functions \code{betapriorGauss}, \code{omegapriorGauss}, \code{etapriorGauss}, 
##' \code{indepGaussianprior} and \code{derivindepGaussianprior}: the first three of which provide a mechanism for storing and retrieving the 
##' parameters of the priors; the fourth, a function for evaluating the log of the prior for a given set of parameter values; and the fifth, a 
##' function for evaluating the first and second derivatives of the log of the prior. It is assumed that parameters are a priori independent. 
##' The user interested in using other priors is encouraged to look at the structure of the five functions mentioned above.
##' 
##' @param betaprior prior for beta, the covariate effects
##' @param omegaprior prior for omega, the parameters of the baseline hazard
##' @param etaprior prior for eta, the parameters of the latent field
##' @param call function to evaluate the log-prior e.g. logindepGaussianprior
##' @param derivative function to evaluate the first and second derivatives of the prior 
##' @return an object of class mcmcPriors
##' @seealso \link{survspat}, \link{betapriorGauss}, \link{omegapriorGauss}, \link{etapriorGauss}, \link{indepGaussianprior}, \link{derivindepGaussianprior}
##' @export

mcmcPriors <- function(betaprior=NULL,omegaprior=NULL,etaprior=NULL,call=NULL,derivative=NULL){
    retlist <- list()
    retlist$betaprior <- betaprior
    retlist$omegaprior <- omegaprior
    retlist$etaprior <- etaprior
    retlist$call <- call
    retlist$derivative <- derivative
    class(retlist) <- "mcmcPriors"
    return(retlist)
}



##' betapriorGauss function
##'
##' A function to define Gaussian priors for beta. This function simply stores a vector of means and standard deviations to be passed to the main MCMC function, survspat.
##'
##' @param mean the prior mean, a vector of length 1 or more. 1 implies a common mean.
##' @param sd the prior standard deviation, a vector of length 1 or more. 1 implies a common standard deviation.
##' @return an object of class "betapriorGauss"
##' @seealso \link{survspat}, \link{betapriorGauss}, \link{omegapriorGauss}, \link{etapriorGauss}, \link{indepGaussianprior}, \link{derivindepGaussianprior}
##' @export

betapriorGauss <- function(mean,sd){
    retlist <- list()
    retlist$mean <- mean
    retlist$sd <- sd
    class(retlist) <- "betapriorGauss"
    return(retlist) 
}



##' omegapriorGauss function
##'
##' A function to define Gaussian priors for omega. This function simply stores a vector of means and standard deviations to be passed to the main MCMC function, survspat.
##'
##' @param mean the prior mean, a vector of length 1 or more. 1 implies a common mean.
##' @param sd the prior standard deviation, a vector of length 1 or more. 1 implies a common standard deviation.
##' @return an object of class "omegapriorGauss"
##' @seealso \link{survspat}, \link{betapriorGauss}, \link{omegapriorGauss}, \link{etapriorGauss}, \link{indepGaussianprior}, \link{derivindepGaussianprior}
##' @export

omegapriorGauss <- function(mean,sd){
    retlist <- list()
    retlist$mean <- mean
    retlist$sd <- sd
    class(retlist) <- "omegapriorGauss"
    return(retlist) 
}



##' etapriorGauss function
##'
##' A function to define Gaussian priors for eta. This function simply stores a vector of means and standard deviations to be passed to the main MCMC function, survspat.
##'
##' @param mean the prior mean, a vector of length 1 or more. 1 implies a common mean.
##' @param sd the prior standard deviation, a vector of length 1 or more. 1 implies a common standard deviation.
##' @return an object of class "etapriorGauss"
##' @seealso \link{survspat}, \link{betapriorGauss}, \link{omegapriorGauss}, \link{etapriorGauss}, \link{indepGaussianprior}, \link{derivindepGaussianprior}
##' @export

etapriorGauss <- function(mean,sd){
    retlist <- list()
    retlist$mean <- mean
    retlist$sd <- sd
    class(retlist) <- "etapriorGauss"
    return(retlist) 
}




##' indepGaussianprior function
##'
##' A function for evaluating the log of an independent Gaussian prior for a given set of parameter values.
##'
##' @param beta parameter beta at which prior is to be evaluated 
##' @param omega parameter omega at which prior is to be evaluated
##' @param eta parameter eta at which prior is to be evaluated
##' @param priors an object of class mcmcPriors, see ?mcmcPriors
##' @return the log of the prior evaluated at the given parameter values
##' @seealso \link{survspat}, \link{betapriorGauss}, \link{omegapriorGauss}, \link{etapriorGauss}, \link{indepGaussianprior}, \link{derivindepGaussianprior}
##' @export

indepGaussianprior <- function(beta=NULL,omega=NULL,eta=NULL,priors){
    
    lp <- 0
    if(!is.null(priors$betaprior)){
        lp <- lp + sum(dnorm(beta,priors$betaprior$mean,priors$betaprior$sd,log=TRUE))
    }

    if(!is.null(priors$omegaprior)){
        lp <- lp + sum(dnorm(omega,priors$omegaprior$mean,priors$omegaprior$sd,log=TRUE))
    }
    
    if(!is.null(priors$etaprior)){
        lp <- lp + sum(dnorm(eta,priors$etaprior$mean,priors$etaprior$sd,log=TRUE))
    }
    
    return(lp)
}


##' derivindepGaussianprior function
##'
##' A function for evaluating the first and second derivatives of the log of an independent Gaussian prior
##'
##' @param beta a vector, the parameter beta
##' @param omega a vector, the parameter omega
##' @param eta a vector, the parameter eta 
##' @param priors an object of class 'mcmcPrior', see ?mcmcPrior
##' @return returns the first and second derivatives of the prior
##' @seealso \link{survspat}, \link{betapriorGauss}, \link{omegapriorGauss}, \link{etapriorGauss}, \link{indepGaussianprior}, \link{derivindepGaussianprior}
##' @export

derivindepGaussianprior <- function(beta=NULL,omega=NULL,eta=NULL,priors){
    if(is.null(beta)){
        deriv1beta <- c()
        deriv2beta <- c()
    }
    else{
        deriv1beta <- (-1/priors$betaprior$sd^2)*(beta-priors$betaprior$mean)
        sdbeta <- priors$betaprior$sd
        if (length(priors$betaprior$sd)<length(beta)){
            sdbeta <- rep(priors$betaprior$sd,length(beta))
        }
        deriv2beta <- -1/sdbeta^2
    }
    if(is.null(omega)){
        deriv1omega <- c()
        deriv2omega <- c()
    }
    else{
        deriv1omega <- (-1/priors$omegaprior$sd^2)*(omega-priors$omegaprior$mean)
        sdomega <- priors$omegaprior$sd 
        if (length(priors$omegaprior$sd)<length(omega)){
            sdomega <- rep(priors$omegaprior$sd,length(omega))
        } 
        deriv2omega <- -1/sdomega^2 
    }
    if(is.null(eta)){
        deriv1eta <- c()
        deriv2eta <- c()
    }
    else{
        deriv1eta <- (-1/priors$etaprior$sd^2)*(eta-priors$etaprior$mean)
        sdeta <- priors$etaprior$sd
        if (length(priors$etaprior$sd)<length(eta)){
            sdeta <- rep(priors$etaprior$sd,length(eta))
        }
        deriv2eta <- -1/sdeta^2
    }

    deriv1 <- c(deriv1beta,deriv1omega,deriv1eta)
    deriv2 <- diag(c(deriv2beta,deriv2omega,deriv2eta)) # in fact, the 2nd derivative with respect to eta is not necessary, as this is dealt with elsewhere.
    
    return(list(deriv1=deriv1,deriv2=deriv2))
}