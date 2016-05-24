##' distinfo function
##'
##' Generic function for returning information about the class of baseline hazard functions employed.
##'
##' @param obj an object
##' @param ... additional argument -- currently there are none, but this is for extensibility
##' @return method distinfo
##' @seealso \link{distinfo.basehazardspec}, \link{exponentialHaz}, \link{weibullHaz}, \link{gompertzHaz}, \link{makehamHaz}, \link{tpowHaz} 
##' @export

distinfo <- function(obj,...){
    UseMethod("distinfo")
}



##' basehazard function
##'
##' Generic function for computing the baseline hazard
##'
##' @param obj an object
##' @param ... additional arguments -- currently there are none, but this is for extensibility
##' @return method basehazard
##' @seealso \link{basehazard.basehazardspec}, \link{exponentialHaz}, \link{weibullHaz}, \link{gompertzHaz}, \link{makehamHaz}, \link{tpowHaz} 
##' @export

basehazard <- function(obj,...){
    UseMethod("basehazard")
}



##' gradbasehazard function
##'
##' Generic function for computing the gradient of the baseline hazard
##'
##' @param obj an object
##' @param ... additional arguments -- currently there are none, but this is for extensibility
##' @return method gradbasehazard
##' @seealso \link{gradbasehazard.basehazardspec}, \link{exponentialHaz}, \link{weibullHaz}, \link{gompertzHaz}, \link{makehamHaz}, \link{tpowHaz} 
##' @export

gradbasehazard <- function(obj,...){
    UseMethod("gradbasehazard")
}



##' hessbasehazard function
##'
##' Generic function for computing the hessian of the baseline hazard
##'
##' @param obj an object
##' @param ... additional arguments -- currently there are none, but this is for extensibility
##' @return method hessbasehazard
##' @seealso \link{hessbasehazard.basehazardspec}, \link{exponentialHaz}, \link{weibullHaz}, \link{gompertzHaz}, \link{makehamHaz}, \link{tpowHaz} 
##' @export

hessbasehazard <- function(obj,...){
    UseMethod("hessbasehazard")
}



##' cumbasehazard function
##'
##' Generic function for computing the cumulative baseline hazard
##'
##' @param obj an object
##' @param ... additional arguments -- currently there are none, but this is for extensibility
##' @return method cumbasehazard
##' @seealso \link{cumbasehazard.basehazardspec}, \link{exponentialHaz}, \link{weibullHaz}, \link{gompertzHaz}, \link{makehamHaz}, \link{tpowHaz} 
##' @export

cumbasehazard <- function(obj,...){
    UseMethod("cumbasehazard")
}



##' gradcumbasehazard function
##'
##' Generic function for computing the gradient of the cumulative baseline hazard 
##'
##' @param obj an object
##' @param ... additional arguments -- currently there are none, but this is for extensibility
##' @return method gradcumbasehazard
##' @seealso \link{gradcumbasehazard.basehazardspec}, \link{exponentialHaz}, \link{weibullHaz}, \link{gompertzHaz}, \link{makehamHaz}, \link{tpowHaz} 
##' @export

gradcumbasehazard <- function(obj,...){
    UseMethod("gradcumbasehazard")
}



##' hesscumbasehazard function
##'
##' Generic function for computing the Hessian of the cumulative baseline hazard 
##'
##' @param obj an object
##' @param ... additional arguments -- currently there are none, but this is for extensibility
##' @return method hesscumbasehazard
##' @seealso \link{hesscumbasehazard.basehazardspec}, \link{exponentialHaz}, \link{weibullHaz}, \link{gompertzHaz}, \link{makehamHaz}, \link{tpowHaz} 
##' @export

hesscumbasehazard <- function(obj,...){
    UseMethod("hesscumbasehazard")
}



##' densityquantile function
##'
##' Generic function for computing quantiles of the density function for a given baseline hazard. This may not be analytically tractable.
##'
##' @param obj an object
##' @param ... additional arguments -- currently there are none, but this is for extensibility
##' @return method densityquantile
##' @seealso \link{densityquantile.basehazardspec}, \link{exponentialHaz}, \link{weibullHaz}, \link{gompertzHaz}, \link{makehamHaz}, \link{tpowHaz}
##' @export

densityquantile <- function(obj,...){
    UseMethod("densityquantile")
}



##' distinfo.basehazardspec function
##'
##' A function to retrive information on the baseline hazard distribution of choice 
##'
##' @param obj an object of class basehazardspec
##' @param ... additional arguments -- currently there are none, but this is for extensibility 
##' @return a function returning information on the baseline hazard distribution of choice
##' @seealso \link{exponentialHaz}, \link{weibullHaz}, \link{gompertzHaz}, \link{makehamHaz}, \link{tpowHaz}
##' @export

distinfo.basehazardspec <- function(obj,...){
    return(obj$distinfo)
}



##' basehazard.basehazardspec function
##'
##' A function to retrieve the baseline hazard function
##'
##' @param obj an object of class basehazardspec
##' @param ... additional arguments -- currently there are none, but this is for extensibility
##' @return a function returning the baseline hazard
##' @seealso \link{exponentialHaz}, \link{weibullHaz}, \link{gompertzHaz}, \link{makehamHaz}, \link{tpowHaz}
##' @export

basehazard.basehazardspec <- function(obj,...){
    return(obj$basehazard)
}



##' gradbasehazard.basehazardspec function
##'
##' A function to retrieve the gradient of the baseline hazard function
##'
##' @param obj an object of class basehazardspec 
##' @param ... additional arguments -- currently there are none, but this is for extensibility
##' @return a function returning the gradient of the baseline hazard
##' @seealso \link{exponentialHaz}, \link{weibullHaz}, \link{gompertzHaz}, \link{makehamHaz}, \link{tpowHaz}
##' @export

gradbasehazard.basehazardspec <- function(obj,...){
    return(obj$gradbasehazard)
}



##' hessbasehazard.basehazardspec function
##'
##' A function to retrieve the Hessian of the baseline hazard function 
##'
##' @param obj an object of class basehazardspec
##' @param ... additional arguments -- currently there are none, but this is for extensibility
##' @return a function returning the Hessian of the baseline hazard
##' @seealso \link{exponentialHaz}, \link{weibullHaz}, \link{gompertzHaz}, \link{makehamHaz}, \link{tpowHaz}
##' @export

hessbasehazard.basehazardspec <- function(obj,...){
    return(obj$hessbasehazard)
}



##' cumbasehazard.basehazardspec function
##'
##' A function to retrieve the cumulative baseline hazard function
##'
##' @param obj an object of class basehazardspec 
##' @param ... additional arguments -- currently there are none, but this is for extensibility
##' @return a function returning the cumulative baseline hazard
##' @seealso \link{exponentialHaz}, \link{weibullHaz}, \link{gompertzHaz}, \link{makehamHaz}, \link{tpowHaz}
##' @export

cumbasehazard.basehazardspec <- function(obj,...){
    return(obj$cumbasehazard)
}



##' gradcumbasehazard.basehazardspec function
##'
##' A function to retrieve the gradient of the cumulative baseline hazard function
##'
##' @param obj an object of class basehazardspec 
##' @param ... additional arguments -- currently there are none, but this is for extensibility
##' @return a function returning the gradient of the cumulative baseline hazard 
##' @seealso \link{exponentialHaz}, \link{weibullHaz}, \link{gompertzHaz}, \link{makehamHaz}, \link{tpowHaz}
##' @export

gradcumbasehazard.basehazardspec <- function(obj,...){
    return(obj$gradcumbasehazard)
}



##' hesscumbasehazard.basehazardspec function
##'
##' A function to retrieve the hessian of the cumulative baseline hazard function
##'
##' @param obj an object of class basehazardspec 
##' @param ... additional arguments -- currently there are none, but this is for extensibility
##' @return a function returning the hessian of the cumulative baseline hazard
##' @seealso \link{exponentialHaz}, \link{weibullHaz}, \link{gompertzHaz}, \link{makehamHaz}, \link{tpowHaz}
##' @export

hesscumbasehazard.basehazardspec <- function(obj,...){
    return(obj$hesscumbasehazard)
}



##' densityquantile.basehazardspec function
##'
##' A function to retrieve the quantiles of the density function
##'
##' @param obj an object of class basehazardspec 
##' @param ... additional arguments -- currently there are none, but this is for extensibility
##' @return a function returning the density quantiles
##' @seealso \link{exponentialHaz}, \link{weibullHaz}, \link{gompertzHaz}, \link{makehamHaz}, \link{tpowHaz}
##' @export

densityquantile.basehazardspec <- function(obj,...){
    return(obj$densityquantile)
}
