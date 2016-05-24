###
# temporalAtRisk class definition and functions
###

##' temporalAtRisk function
##'
##' Generic function used in the construction of temporalAtRisk objects. A temporalAtRisk object describes the at risk
##' population globally in an observation time window [t_1,t_2]. Therefore, for any t in [t_1,t_2], a temporalAtRisk object should be able
##' to return the global at risk population, mu(t) = E(number of cases in the unit time interval containing t). This is in contrast to the class of 
##' \link{spatialAtRisk} objects, which describe the spatial inhomogeneity in the population at risk, lambda(s).
##'
##' Note that in the prediction routine, \link{lgcpPredict}, and the simulation routine, \link{lgcpSim}, time discretisation is achieved
##' using \code{as.integer} on both observation times and time limits t_1 and t_2 (which may be stored as non-integer values). The
##' functions that create temporalAtRisk objects therefore return piecewise cconstant step-functions. that can be evaluated for any real
##' t in [t_1,t_2], but with the restriction that mu(t_i) = mu(t_j) whenever \code{as.integer(t_i)==as.integer(t_j)}.
##'
##' A temporalAtRisk object may be (1) 'assumed known', or (2) scaled to a particular dataset. In the latter case, in the routines available (\link{temporalAtRisk.numeric}
##' and \link{temporalAtRisk.function}), the stppp dataset of interest should be referenced, in which case the scaling of mu(t) will be done 
##' automatically. Otherwise, for example for simulation purposes, no scaling of mu(t) occurs, and it is assumed that the mu(t) corresponds to the
##' expected number of cases during the unit time interval containnig t.
##' For reference purposes, the following is a mathematical description of a log-Gaussian Cox Process, it is best viewed in the pdf version of the manual.
##'
##' Let \eqn{\mathcal Y(s,t)}{\mathcal Y(s,t)} be a spatiotemporal Gaussian process, \eqn{W\subset R^2}{W\subset R^2} be an 
##' observation window in space and \eqn{T\subset R_{\geq 0}}{T\subset R_{\geq 0}} be an interval of time of interest. 
##' Cases occur at spatio-temporal positions \eqn{(x,t) \in W \times T}{(x,t) \in W \times T} 
##'  according to an inhomogeneous spatio-temporal Cox process,
##' i.e. a Poisson process with a stochastic intensity \eqn{R(x,t)}{R(x,t)},
##'   The number of cases, \eqn{X_{S,[t_1,t_2]}}{X_{S,[t_1,t_2]}}, arising in 
##'   any \eqn{S \subseteq W}{S \subseteq W} during the interval \eqn{[t_1,t_2]\subseteq T}{[t_1,t_2]\subseteq T} is 
##'   then Poisson distributed conditional on \eqn{R(\cdot)}{R(\cdot)},
##' \deqn{X_{S,[t_1,t_2]} \sim \mbox{Poisson}\left\{\int_S\int_{t_1}^{t_2} R(s,t)d sd t\right\}}{X_{S,[t_1,t_2]} \sim \mbox{Poisson}\left\{\int_S\int_{t_1}^{t_2} R(s,t)d sd t\right\}.}
##' Following Brix and Diggle (2001) and Diggle et al (2005), the intensity is decomposed multiplicatively as
##' \deqn{R(s,t) = \lambda(s)\mu(t)\exp\{\mathcal Y(s,t)\}.}{R(s,t) = \lambda(s)\mu(t)Exp\{\mathcal Y(s,t)\}.}
##' In the above, the fixed spatial component, \eqn{\lambda:R^2\mapsto R_{\geq 0}}{\lambda:R^2\mapsto R_{\geq 0}}, 
##' is a known function, proportional to the population at risk at each point in space and scaled so that
##' \deqn{\int_W\lambda(s)d s=1,}{\int_W\lambda(s)d s=1,}
##' whilst the fixed temporal component, 
##'  \eqn{\mu:R_{\geq 0}\mapsto R_{\geq 0}}{\mu:R_{\geq 0}\mapsto R_{\geq 0}}, is also a known function with
##' \deqn{\mu(t) \delta t = E[X_{W,\delta t}],}{\mu(t) \delta t = E[X_{W,\delta t}],}
##' for \eqn{t}{t} in a small interval of time, \eqn{\delta t}{\delta t}, over which the rate of the process over \eqn{W}{W} can be considered constant.
##'
##' @param obj an object
##' @param ... additional arguments
##' @return method temporalAtRisk
##' \enumerate{
##'     \item Brix A, Diggle PJ (2001). Spatiotemporal Prediction for log-Gaussian Cox processes. Journal of the Royal Statistical Society, Series B, 63(4), 823-841.
##'     \item Diggle P, Rowlingson B, Su T (2005). Point Process Methodology for On-line Spatio-temporal Disease Surveillance. Environmetrics, 16(5), 423-434.
##' }
##' @seealso \link{spatialAtRisk}, \link{lgcpPredict}, \link{lgcpSim}, \link{temporalAtRisk.numeric}, \link{temporalAtRisk.function}, \link{constantInTime}, \link{constantInTime.numeric}, \link{constantInTime.stppp}, \link{print.temporalAtRisk}, \link{plot.temporalAtRisk}
##' @export

temporalAtRisk <- function(obj,...){
    UseMethod("temporalAtRisk")
}



##' temporalAtRisk.numeric function
##'
##' Create a temporalAtRisk object from a numeric vector.
##'
##' Note that in the prediction routine, \link{lgcpPredict}, and the simulation routine, \link{lgcpSim}, time discretisation is achieved
##' using \code{as.integer} on both observation times and time limits t_1 and t_2 (which may be stored as non-integer values). The
##' functions that create temporalAtRisk objects therefore return piecewise constant step-functions that can be evaluated for any real
##' t in [t_1,t_2], but with the restriction that mu(t_i) = mu(t_j) whenever \code{as.integer(t_i)==as.integer(t_j)}.
##'
##' A temporalAtRisk object may be (1) 'assumed known', corresponding to the default argument \code{xyt=NULL}; or (2) scaled to a particular dataset
##' (argument xyt=[stppp object of interest]). In the latter case, in the routines available (\link{temporalAtRisk.numeric}
##' and \link{temporalAtRisk.function}), the dataset of interest should be referenced, in which case the scaling of mu(t) will be done 
##' automatically. Otherwise, for example for simulation purposes, no scaling of mu(t) occurs, and it is assumed that the mu(t) corresponds to the
##' expected number of cases during the unit time interval containing t.
##'
##' @method temporalAtRisk numeric
##' @param obj a numeric vector of length (tlim[2]-tlim[1] + 1) giving the temporal intensity up to a constant of proportionality at each integer time within the interval defined by tlim
##' @param tlim an integer vector of length 2 giving the time limits of the observation window
##' @param xyt an object of class stppp. If NULL (default) then the function returned is not scaled. Otherwise, the function is scaled so that f(t) = expected number of counts at time t.
##' @param warn Issue a warning if the given temporal intensity treated is treated as 'known'?
##' @param ... additional arguments
##' @return a function f(t) giving the temporal intensity at time t for integer t in the interval as.integer([tlim[1],tlim[2]]) of class temporalAtRisk
##' \enumerate{
##'     \item Brix A, Diggle PJ (2001). Spatiotemporal Prediction for log-Gaussian Cox processes. Journal of the Royal Statistical Society, Series B, 63(4), 823-841.
##'     \item Diggle P, Rowlingson B, Su T (2005). Point Process Methodology for On-line Spatio-temporal Disease Surveillance. Environmetrics, 16(5), 423-434.
##' }
##' @seealso \link{temporalAtRisk}, \link{spatialAtRisk},  \link{temporalAtRisk.function}, \link{constantInTime}, \link{constantInTime.numeric}, \link{constantInTime.stppp}, \link{print.temporalAtRisk}, \link{plot.temporalAtRisk}
##' @export

temporalAtRisk.numeric <- function(obj,tlim,xyt=NULL,warn=TRUE,...){
    tlim <- sort(as.integer(tlim))
    tvec <- tlim[1]:tlim[2]
    if (length(obj)!=length(tvec)){
        stop(paste("obj must be a vector of length",length(tvec)))
    }
    scale <- 1    
    if (!is.null(xyt)){
        if (!all(xyt$tlim==tlim)){
            stop("xyt$tlim and tlim incompatible.")
        }
        scale <- xyt$n / sum(obj) # maximum likelihood estimate of alpha in X_i~Poisson(alpha mu_i), where mu_i are known (and equal to obj in this case)        
    }
    else{
        if(warn){
            warning("No xyt specified, given temporal intensity treated is treated as known")
        }
    }
    ft <- function(t){
        if(length(t)>1){
            stop("Function only works for scalar numeric objects, ie a vector of length 1.")
        }
        if (!any(as.integer(t)==tvec)){
            return(NA)
        }
        return(obj[which(as.integer(t)==tvec)] * scale)
    }
    object <- ft
    attr(object,"tlim") <- tlim
    class(object) <- c("temporalAtRisk","function")
    return(object) 
}



##' temporalAtRisk.function function
##'
##' Create a temporalAtRisk object from a function.
##'
##' Note that in the prediction routine, \link{lgcpPredict}, and the simulation routine, \link{lgcpSim}, time discretisation is achieved
##' using \code{as.integer} on both observation times and time limits t_1 and t_2 (which may be stored as non-integer values). The
##' functions that create temporalAtRisk objects therefore return piecewise cconstant step-functions. that can be evaluated for any real
##' t in [t_1,t_2], but with the restriction that mu(t_i) = mu(t_j) whenever \code{as.integer(t_i)==as.integer(t_j)}.
##'
##' A temporalAtRisk object may be (1) 'assumed known', corresponding to the default argument \code{xyt=NULL}; or (2) scaled to a particular dataset
##' (argument xyt=[stppp object of interest]). In the latter case, in the routines available (\link{temporalAtRisk.numeric}
##' and \link{temporalAtRisk.function}), the dataset of interest should be referenced, in which case the scaling of mu(t) will be done 
##' automatically. Otherwise, for example for simulation purposes, no scaling of mu(t) occurs, and it is assumed that the mu(t) corresponds to the
##' expected number of cases during the unit time interval containnig t.
##'
##' @method temporalAtRisk function
##' @param obj a function accepting single, scalar, numeric argument, t, that returns the temporal intensity for time t
##' @param tlim an integer vector of length 2 giving the time limits of the observation window
##' @param xyt an object of class stppp. If NULL (default) then the function returned is not scaled. Otherwise, the function is scaled so that f(t) = expected number of counts at time t.
##' @param warn Issue a warning if the given temporal intensity treated is treated as 'known'?
##' @param ... additional arguments
##' @return a function f(t) giving the temporal intensity at time t for integer t in the interval [tlim[1],tlim[2]] of class temporalAtRisk
##' \enumerate{
##'     \item Brix A, Diggle PJ (2001). Spatiotemporal Prediction for log-Gaussian Cox processes. Journal of the Royal Statistical Society, Series B, 63(4), 823-841.
##'     \item Diggle P, Rowlingson B, Su T (2005). Point Process Methodology for On-line Spatio-temporal Disease Surveillance. Environmetrics, 16(5), 423-434.
##' }
##' @seealso \link{temporalAtRisk}, \link{spatialAtRisk}, \link{temporalAtRisk.numeric},  \link{constantInTime}, \link{constantInTime.numeric}, \link{constantInTime.stppp}, \link{print.temporalAtRisk}, \link{plot.temporalAtRisk}
##' @export

temporalAtRisk.function <- function(obj,tlim,xyt=NULL,warn=TRUE,...){
    tlim <- sort(as.integer(tlim))
    scale <- 1
    object <- list()
    if (!is.null(xyt)){
        if (!all(xyt$tlim==tlim)){
            stop("xyt$tlim and tlim incompatible.")
        }
        scale <- xyt$n / sum(sapply(tlim[1]:tlim[2],obj)) # maximum likelihood estimate of alpha in X_i~Poisson(alpha mu_i), where mu_i are known (and equal to obj in this case)        
    }
    else{
        if(warn){
            warning("No xyt specified, given temporal intensity treated is treated as known")
        }
    }
    ft <- function(t){
        if(length(t)>1){
            stop("obj must be a function of a numeric scalar.")
        }
        if (as.integer(t)<tlim[1] | as.integer(t)>tlim[2]){
            return(NA)
        }
        return(obj(as.integer(t)) * scale)
    }
    object <- ft
    attr(object,"tlim") <- tlim
    class(object) <- c("temporalAtRisk","function")
    return(object)
}



##' constantInTime function
##'
##' Generic function for creating constant-in-time temporalAtRisk objects, that is for models where mu(t) can be assumed to
##' be constant in time. The assumption being that the global at-risk population does not change in size over time.
##'
##' For further details of temporalAtRisk objects, see ?temporalAtRisk>
##'
##' @param obj an object
##' @param ... additional arguments
##' @return method constantInTime
##' @seealso \link{temporalAtRisk}, \link{spatialAtRisk}, \link{temporalAtRisk.numeric}, \link{temporalAtRisk.function}, \link{constantInTime.numeric}, \link{constantInTime.stppp}, \link{print.temporalAtRisk}, \link{plot.temporalAtRisk}
##' @export

constantInTime <- function(obj,...){
    UseMethod("constantInTime")
}



##' constantInTime.numeric function
##'
##' Create a constant-in-time temporalAtRisk object from a numeric object of length 1. The returned temporalAtRisk object is assumed to have
##' been scaled correctly by the user so that mu(t) = E(number of cases in a unit time interval).
##'
##' For further details of temporalAtRisk objects, see ?temporalAtRisk>
##'
##' @method constantInTime numeric
##' @param obj numeric constant
##' @param tlim vector of length 2 giving time limits
##' @param warn Issue a warning if the given temporal intensity treated is treated as 'known'?
##' @param ... additional arguments
##' @return a function f(t) giving the (constant) temporal intensity at time t for integer t in the interval [tlim[1],tlim[2]] of class temporalAtRisk
##' @seealso \link{temporalAtRisk}, \link{spatialAtRisk}, \link{temporalAtRisk.numeric}, \link{temporalAtRisk.function}, \link{constantInTime}, \link{constantInTime.stppp}, \link{print.temporalAtRisk}, \link{plot.temporalAtRisk}, 
##' @export

constantInTime.numeric <- function(obj,tlim,warn=TRUE,...){
    if (length(obj)>1){
        stop("obj must be a single real number")
    }
    f <- function(t){
        return(obj)
    }
    return(temporalAtRisk(f,tlim,warn=warn))
}



##' constantInTime.stppp function
##'
##' Create a constant-in-time temporalAtRisk object from an stppp object. The returned temporalAtRisk object is scaled to return
##' mu(t) = E(number of cases in a unit time interval).
##'
##' For further details of temporalAtRisk objects, see ?temporalAtRisk>
##'
##' @method constantInTime stppp
##' @param obj an object of class stppp.
##' @param ... additional arguments
##' @return a function f(t) giving the (constant) temporal intensity at time t for integer t in the interval [tlim[1],tlim[2]] of class temporalAtRisk
##' @seealso \link{temporalAtRisk}, \link{spatialAtRisk}, \link{temporalAtRisk.numeric}, \link{temporalAtRisk.function}, \link{constantInTime}, \link{constantInTime.numeric}, \link{print.temporalAtRisk}, \link{plot.temporalAtRisk}, 
##' @export

constantInTime.stppp <- function(obj,...){
    f <- function(t){
        return(1)
    }
    return(temporalAtRisk(obj=f,tlim=obj$tlim,xyt=obj))
}



##' print.temporalAtRisk function
##'
##' Printing method for temporalAtRisk objects.
##'
##' @method print temporalAtRisk
##' @param x an object
##' @param ... additional arguments
##' @return print the object
##' @seealso \link{temporalAtRisk}, \link{spatialAtRisk}, \link{temporalAtRisk.numeric}, \link{temporalAtRisk.function}, \link{constantInTime}, \link{constantInTime.numeric}, \link{constantInTime.stppp}, \link{plot.temporalAtRisk}
##' @export

print.temporalAtRisk <- function(x,...){
    cat("temporalAtRisk object\n")
    NextMethod("print",x)
    cat(paste("   Time Window : [",attr(x,"tlim")[1],",",attr(x,"tlim")[2],"]\n"))
}



##' plot.temporalAtRisk function
##'
##' Pot a temporalAtRisk object.
##'
##' @method plot temporalAtRisk
##' @param x an object
##' @param ... additional arguments
##' @return print the object
##' @seealso \link{temporalAtRisk}, \link{spatialAtRisk}, \link{temporalAtRisk.numeric}, \link{temporalAtRisk.function}, \link{constantInTime}, \link{constantInTime.numeric}, \link{constantInTime.stppp}, \link{print.temporalAtRisk}, 
##' @export

plot.temporalAtRisk <- function(x,...){
    tlim <- attr(x,"tlim")
    time <- seq(tlim[1],tlim[2],length.out=10000)
    at_risk <- sapply(time,x)
    plot(time,at_risk,type="l",...)
}
