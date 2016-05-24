###
# Generic functions for the class adaptivemcmc    
###



##' initialiseAMCMC function
##'
##' A generic to be used for the purpose of user-defined adaptive MCMC schemes,
##' initialiseAMCMC tells the MALA algorithm which value of h to use first. See
##' lgcp vignette, code{vignette("lgcp")}, for further details on writing adaptive MCMC schemes.
##'
##' @param obj an object    
##' @param ... additional arguments  
##' @return method intialiseAMCMC
##' @seealso \link{initialiseAMCMC.constanth}, \link{initialiseAMCMC.andrieuthomsh}
##' @export

initialiseAMCMC <- function(obj,...){
    UseMethod("initialiseAMCMC")
}



##' updateAMCMC function
##'
##' A generic to be used for the purpose of user-defined adaptive MCMC schemes,
##' updateAMCMC tells the MALA algorithm how to update the value of h. See
##' lgcp vignette, code{vignette("lgcp")}, for further details on writing adaptive MCMC schemes.
##'
##' @param obj an object    
##' @param ... additional arguments  
##' @return method updateAMCMC
##' @seealso \link{updateAMCMC.constanth}, \link{updateAMCMC.andrieuthomsh}
##' @export

updateAMCMC <- function(obj,...){
    UseMethod("updateAMCMC")
}



###
# Base functions for creating adaptive schemes
###



##' constanth function
##'
##' This function is used to set up a constant acceptance scheme in the argument 
##' \code{mcmc.control} of the function \code{lgcpPredict}. The scheme is set via
##' the \code{mcmcpars} function.
##'
##' @param h an object 
##' @return object of class constanth
##' @seealso \link{mcmcpars}, \link{lgcpPredict}
##' @examples constanth(0.01)
##' @export

constanth <- function(h){
    obj <- h
    class(obj) <- c("constanth","adaptivemcmc")
    return(obj)
}



##' andrieuthomsh function
##'
##' A Robbins-Munro stochastic approximation update is used to adapt the tuning parameter of the proposal kernel. 
##' The idea is to update the tuning parameter at each iteration of the sampler:
##' \deqn{h^{(i+1)} = h^{(i)} + \eta^{(i+1)}(\alpha^{(i)} - \alpha_{opt}),}{h^{(i+1)} = h^{(i)} + \eta^{(i+1)}(\alpha^{(i)} - \alpha_{opt}),}
##' where \eqn{h^{(i)}}{h^{(i)}} and \eqn{\alpha^{(i)}}{\alpha^{(i)}} are the tuning parameter and acceptance probability at iteration 
##' \eqn{i}{i} and \eqn{\alpha_{opt}}{\alpha_{opt}} is a target acceptance probability. For Gaussian targets, and in the limit 
##' as the dimension of the problem tends to infinity, an appropriate target acceptance probability for 
##' MALA algorithms is 0.574. The sequence \eqn{\{\eta^{(i)}\}}{\{\eta^{(i)}\}} is chosen so that 
##' \eqn{\sum_{i=0}^\infty\eta^{(i)}}{\sum_{i=0}^\infty\eta^{(i)}} is infinite whilst \eqn{\sum_{i=0}^\infty\left(\eta^{(i)}\right)^{1+\epsilon}}{\sum_{i=0}^\infty\left(\eta^{(i)}\right)^{1+\epsilon}} is 
##' finite for \eqn{\epsilon>0}{\epsilon>0}. These two conditions ensure that any value of \eqn{h}{h} can be reached, but in a way that 
##' maintains the ergodic behaviour of the chain. One class of sequences with this property is,
##' \deqn{\eta^{(i)} = \frac{C}{i^\alpha},}{\eta^{(i)} = \frac{C}{i^\alpha},}
##' where \eqn{\alpha\in(0,1]}{\alpha\in(0,1]} and \eqn{C>0}{C>0}.The scheme is set via
##' the \code{mcmcpars} function.
##'
##' @param inith initial h
##' @param alpha parameter \eqn{\alpha}{\alpha}
##' @param C parameter \eqn{C}{C}
##' @param targetacceptance target acceptance probability
##' @return an object of class andrieuthomsh
##' @references \enumerate{
##'     \item Andrieu C, Thoms J (2008). A tutorial on adaptive MCMC. Statistics and Computing, 18(4), 343-373.
##'     \item Robbins H, Munro S (1951). A Stochastic Approximation Methods. The Annals of Mathematical Statistics, 22(3), 400-407.
##'     \item Roberts G, Rosenthal J (2001). Optimal Scaling for Various Metropolis-Hastings Algorithms. Statistical Science, 16(4), 351-367.
##' }
##' @seealso \link{mcmcpars}, \link{lgcpPredict}
##' @examples andrieuthomsh(inith=1,alpha=0.5,C=1,targetacceptance=0.574)
##' @export

andrieuthomsh <- function(inith,alpha,C,targetacceptance=0.574){
    if (alpha<=0 | alpha>1){
        stop("parameter alpha must be in [0,1]")
    }
    if (C<=0){
        stop("parameter C must be positive")
    }
    obj <- list()
    obj$inith <- inith
    obj$alpha <- alpha
    obj$C <- C
    obj$targetacceptance <- targetacceptance
    
    itno <- 0 # iteration number, gets reset after burnin
    incrit <- function(){
        itno <<- itno + 1
    }
    restit <- function(){
        itno <<- 0
    }
    obj$incritno <- incrit
    obj$restartit <- restit
    
    curh <- inith # at iteration 1, the current value of h is just the initial value
    hupdate <- function(){ 
        curh <<- exp(log(curh) + (C/(itno^alpha))*(get("ac",envir=parent.frame(2))-targetacceptance))
    }
    reth <- function(){
        return(curh)
    }
    obj$updateh <- hupdate
    obj$returncurh <- reth
    
    class(obj) <- c("andrieuthomsh","adaptivemcmc")
    return(obj)
}



###
# Initialisation/updating functions to implement adaptive schemes
###

##' initaliseAMCMC.constanth function
##'
##' Initialises the \link{constanth} adaptive scheme.
##'
##' @method initialiseAMCMC constanth
##' @param obj an object 
##' @param ... additional arguments 
##' @return initial h for scheme
##' @seealso \link{constanth}
##' @export

initialiseAMCMC.constanth <- function(obj,...){
    return(obj)
}



##' initaliseAMCMC.andrieuthomsh function
##'
##' Initialises the \link{andrieuthomsh} adaptive scheme.
##'
##' @method initialiseAMCMC andrieuthomsh
##' @param obj an object 
##' @param ... additional arguments
##' @return initial h for scheme
##' @references \enumerate{
##'     \item Andrieu C, Thoms J (2008). A tutorial on adaptive MCMC. Statistics and Computing, 18(4), 343-373.
##'     \item Robbins H, Munro S (1951). A Stochastic Approximation Methods. The Annals of Mathematical Statistics, 22(3), 400-407.
##'     \item Roberts G, Rosenthal J (2001). Optimal Scaling for Various Metropolis-Hastings Algorithms. Statistical Science, 16(4), 351-367.
##' }
##' @seealso \link{andrieuthomsh}
##' @export

initialiseAMCMC.andrieuthomsh <- function(obj,...){
    return(obj$inith)
}



##' updateAMCMC.constanth function
##'
##' Updates the \link{constanth} adaptive scheme.
##'
##' @method updateAMCMC constanth
##' @param obj an object 
##' @param ... additional arguments
##' @return update and return current h for scheme
##' @seealso \link{constanth}
##' @export

updateAMCMC.constanth <- function(obj,...){
    return(obj)
}



##' updateAMCMC.andrieuthomsh function
##'
##' Updates the \link{andrieuthomsh} adaptive scheme.
##'
##' @method updateAMCMC andrieuthomsh
##' @param obj an object 
##' @param ... additional arguments
##' @return update and return current h for scheme
##' @references \enumerate{
##'     \item Andrieu C, Thoms J (2008). A tutorial on adaptive MCMC. Statistics and Computing, 18(4), 343-373.
##'     \item Robbins H, Munro S (1951). A Stochastic Approximation Methods. The Annals of Mathematical Statistics, 22(3), 400-407.
##'     \item Roberts G, Rosenthal J (2001). Optimal Scaling for Various Metropolis-Hastings Algorithms. Statistical Science, 16(4), 351-367.
##' }
##' @seealso \link{andrieuthomsh}
##' @export

updateAMCMC.andrieuthomsh <- function(obj,...){
    
    mLoop <- get("mcmcloop",envir=parent.frame())
    if(iteration(mLoop)==(mLoop$burnin)+1){
        obj$restartit() # adaptation of h restarted after burnin, this allows h to adapt more quickly to the information learned during burnin
    }    
    obj$incritno() # this line must appear below the above four lines 
        
    obj$updateh()
    return(obj$returncurh())
}

