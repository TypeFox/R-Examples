##' weibullHaz function
##'
##' A function to define a parametric proportional hazards model where the baseline hazard is taken from the Weibull model. 
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
##' @seealso \link{tpowHaz}, \link{exponentialHaz}, \link{gompertzHaz}, \link{makehamHaz}
##' @export

weibullHaz <- function(){

    flist <- list()

    flist$distinfo <- function(){
        retlist <- list()
        retlist$npars <- 2
        retlist$parnames <- c("alpha","lambda")
        retlist$trans <- log
        retlist$itrans <- exp
        retlist$jacobian <- exp
        retlist$hessian <- list(exp,exp)
        return(retlist)
    }
    
    flist$basehazard <- function(pars){
        fun <- function(t){
            return(pars[2]*pars[1]*t^(pars[1]-1)) # in this case alpha=pars[1], lambda=pars[2] 
        }
        return(fun)  
    }
    
    flist$gradbasehazard <- function(pars){
        fun <- function(t){
            return(t^(pars[1]-1)*cbind(pars[2]*(1+pars[1]*log(t)),pars[1])) # in this case alpha=pars[1], lambda=pars[2]
        }
        return(fun)
        
    }
    
    flist$hessbasehazard <- function(pars){
        funfun <- function(t,pars){
            m <- matrix(0,2,2) # note m[2,2]=0 i.e. d2h_0/dlambda^2 = 0
            m[1,2] <- m[2,1] <- t^(pars[1]-1)*(1+pars[1]*log(t))
            m[1,1] <- pars[2]*t^(pars[1]-1)*log(t)*(2+pars[1]*log(t))
            return(m) # in this case alpha=pars[1], lambda=pars[2]
        }
        
        fun <- function(t){
            return(lapply(t,funfun,pars=pars))
        }
        return(fun)
        
    }
    
    flist$cumbasehazard <- function(pars){
        fun <- function(t){
            return(pars[2]*t^(pars[1])) # in this case alpha=pars[1], lambda=pars[2]
        }
        return(fun) 
    }
    
    flist$gradcumbasehazard <- function(pars){
        fun <- function(t){
            return(t^(pars[1])*cbind(pars[2]*log(t),1)) # in this case alpha=pars[1], lambda=pars[2]
        }
        return(fun)    
    }
    
    flist$hesscumbasehazard <- function(pars){
        funfun <- function(t,pars){
            m <- matrix(0,2,2) # note m[2,2]=0 i.e. d2H_0/dlambda^2 = 0
            other <- log(t)*t^pars[1]
            m[1,2] <- m[2,1] <- other 
            m[1,1] <- pars[2]*other*log(t)
            return(m) # in this case alpha=pars[1], lambda=pars[2]
        }
        
        fun <- function(t){
            return(lapply(t,funfun,pars=pars))
        }
        return(fun)
    }
    
    flist$densityquantile <- function(pars,other){
        fun <- function(probs){
            return((-log(1-probs)/(pars[2]*other$expXbetaplusY))^(1/pars[1]))
        }
        return(fun)    
    }


    class(flist) <- c("basehazardspec","list")
    return(flist)
}