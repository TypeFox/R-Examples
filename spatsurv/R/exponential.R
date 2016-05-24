##' exponentialHaz function
##'
##' A function to define a parametric proportional hazards model where the baseline hazard is taken from the exponential model. 
##' This function returns an object inheriting class 'basehazardspec', list of functions 'distinfo', 'basehazard', 'gradbasehazard', 'hessbasehazard',
##' 'cumbasehazard', 'gradcumbasehazard', 'hesscumbasehazard' and 'densityquantile'
##' 
##' The \code{distinfo} function is used to provide basic distribution specific information to other \code{spatsurv} functions. The user is required 
##' to provide the following information in the returned list: \code{npars}, the number of parameters in this distribution; \code{parnames}, 
##' the names of the parameters; \code{trans}, the transformation scale on which the priors will be provided; \code{itrans}, the inverse 
##' transformation function that will be applied to the parameters before the hazard, and other functions are evaluated; \code{jacobian}, 
##' the derivative of the inverse transformation function with respect to each of the parameters; and \code{hessian}, the second derivatives 
##' of the inverse transformation function with respect to each of the parameters -- note that currently the package \code{spatsurv} 
##' only allows the use of functions where the parameters are transformed independently.
##' 
##' The \code{basehazard} function is used to evaluate the baseline hazard function for the distribution of interest. It returns a 
##' function that accepts as input a vector of times, \code{t} and returns a vector.
##' 
##' The \code{gradbasehazard} function is used to evaluate the gradient of the baseline hazard function with respect to the parameters, 
##' this typically returns a vector. It returns a function that accepts as input a vector of times, \code{t}, and returns a matrix.
##' 
##' The \code{hessbasehazard} function is used to evaluate the Hessian of the baseline hazard function. It returns a function that accepts 
##' as input a vector of times, \code{t} and returns a list of hessian matrices corresponding to each \code{t}.
##' 
##' The \code{cumbasehazard} function is used to evaluate the cumulative baseline hazard function for the distribution of interest. 
##' It returns a function that accepts as input a vector of times, \code{t} and returns a vector.
##' 
##' The \code{gradcumbasehazard} function is used to evaluate the gradient of the cumulative baseline hazard function with respect 
##' to the parameters, this typically returns a vector. It returns a function that accepts as input a vector of times, \code{t}, and returns a matrix.
##' 
##' The \code{hesscumbasehazard} function is used to evaluate the Hessian of the cumulative baseline hazard function. It returns a 
##' function that accepts as input a vector of times, \code{t} and returns a list of hessian matrices corresponding to each \code{t}.
##' 
##' The \code{densityquantile} function is used to return quantiles of the density function. This is NOT REQUIRED for running the MCMC, 
##' merely for us in post-processing with the \code{predict} function where \code{type} is 'densityquantile'. In the case of the Weibull 
##' model for the baseline hazard, it can be shown that the q-th quantile is: 
##'
##' @return an object inheriting class 'basehazardspec'
##' @seealso \link{tpowHaz}, \link{gompertzHaz}, \link{makehamHaz}, \link{weibullHaz}
##' @export


exponentialHaz <- function(){

    flist <- list()

    flist$distinfo <- function(){
        retlist <- list()
        retlist$npars <- 1
        retlist$parnames <- "rate"
        retlist$trans <- log
        retlist$itrans <- exp
        retlist$jacobian <- exp
        retlist$hessian <- list(exp)
        return(retlist)
    }
 
    flist$basehazard <- function(pars){
        fun <- function(t){
            return(rep(pars,length(t))) # in this case pars is a 1-vector, the rate    
        }
        return(fun)
    }
 
    flist$gradbasehazard <- function(pars){
        fun <- function(t){
            return(rep(1,length(t))) # in this case pars is a 1-vector, the rate
        }
        return(fun)
    }
   
    flist$hessbasehazard <- function(pars){
        fun <- function(t){
            return(as.list(rep(0,length(t)))) # in this case pars is a 1-vector, the rate
        }
        return(fun)
    }
   
    flist$cumbasehazard <- function(pars){
        fun <- function(t){
            return(pars*t) # in this case pars is a 1-vector, the rate
        }
        return(fun)  
    }
  
    flist$gradcumbasehazard <- function(pars){
        fun <- function(t){
            return(t) # in this case pars is a 1-vector, the rate
        }
        return(fun)    
    }   
    
    flist$hesscumbasehazard <- function(pars){
        fun <- function(t){
            return(as.list(rep(0,length(t)))) # in this case pars is a 1-vector, the rate
        }
        return(fun)
    }  
    
    flist$densityquantile <- function(pars,other){
        fun <- function(probs){
            return(-log(1-probs)/(pars*other$expXbetaplusY))
        }
        return(fun)    
    }

    class(flist) <- c("basehazardspec","list")
    return(flist)
}