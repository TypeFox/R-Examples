##' BsplineHaz function
##'
##' A function to define a parametric proportional hazards model where the baseline hazard is modelled by a basis spline. 
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
##' @param times vector of survival times (both censored and uncensored)
##' @param knots vector of knots in ascending order, must include minimum and maximum values of 'times'
##' @param degree degree of the spline basis, default is 3
##' @param MLinits optional starting values for the non-spatial maximisation routine using optim. Note that we are working with the log of the parameters. Default is -10 for each parameter.
##' @return an object inheriting class 'basehazardspec'
##' @seealso \link{exponentialHaz}, \link{gompertzHaz}, \link{makehamHaz}, \link{weibullHaz} 
##' @export

BsplineHaz <- function(times,knots=quantile(times),degree=3,MLinits=NULL){

	basis <- getBbasis(x=times,knots,degree=degree)

	basismatrix <- Bspline.construct(x=times,basis=basis) 
	basismatrix[basismatrix<0] <- 0 # these are small negative numbers

	cbs <- cumulativeBspline.construct(x=times,basis=basis) 
	
	np <- length(basis$poly)

    flist <- list()

    flist$distinfo <- function(){
        retlist <- list()
        retlist$npars <- np
        retlist$parnames <- paste("lambda",1:np,sep="")
        retlist$trans <- log
        retlist$itrans <- exp
        retlist$jacobian <- exp
        retlist$hessian <- lapply(1:np,function(zz){return(exp)})
        if(is.null(MLinits)){
        	retlist$MLinits <- rep(-10,np)	
        }
        else{
        	retlist$MLinits <- MLinits
        }        
        return(retlist)
    }

    test <- flist$distinfo()
    cat("Using B-spline with ",test$npars," parameters.\n")
    
    flist$basehazard <- function(pars){
        fun <- function(t){
        	idx <- match(t,times)
        	if(any(is.na(idx))){
        		basismatrix <- Bspline.construct(x=t,basis=basis) 
        		basismatrix[basismatrix<0] <- 0 # these are small negative numbers
        		return(colSums(pars*t(basismatrix)))
        	}
        	else{
        		return(colSums(pars*t(basismatrix[idx,,drop=FALSE])))
        	}
        }
        return(fun)  
    }
    
    flist$gradbasehazard <- function(pars){
        fun <- function(t){
        	idx <- match(t,times)
        	return(basismatrix[idx,])
        }
        return(fun)        
    }
    
    flist$hessbasehazard <- function(pars){
        funfun <- function(t,pars){ 
            return(matrix(0,np,np)) 
        }
        
        fun <- function(t){
            return(lapply(t,funfun,pars=pars))
        }
        return(fun)
        
    }
    
    flist$cumbasehazard <- function(pars){
        fun <- function(t){
        	idx <- match(t,times)       	
        	if(any(is.na(idx))){
        		cbs <- cumulativeBspline.construct(x=t,basis=basis) 
        		return(colSums(pars*t(cbs$integral)) + cbs$toadd(pars))
        	}
        	else{
        		return(colSums(pars*t(cbs$integral[idx,,drop=FALSE])) + cbs$toadd(pars))
        	}
        }
        return(fun) 
    }
    
    flist$gradcumbasehazard <- function(pars){
        fun <- function(t){
            idx <- match(t,times)
        	return(cbs$integral[idx,])        
        }
        return(fun)    
    }
    
    flist$hesscumbasehazard <- function(pars){
        funfun <- function(t,pars){
            return(matrix(0,np,np)) 
        }
        
        fun <- function(t){
            return(lapply(t,funfun,pars=pars))
        }
        return(fun)
    }
    
    flist$densityquantile <- function(pars,other){
        fun <- function(probs){
            stop("densityquantile not available yet")
            #return((-log(1-probs)/(pars[2]*other$expXbetaplusY))^(1/pars[1]))
        }
        return(fun)    
    }


    class(flist) <- c("basehazardspec","list")
    return(flist)
}