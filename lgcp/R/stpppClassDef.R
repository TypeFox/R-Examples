###
# stppp class definition and functions
###

##' stppp function
##'
##' Generic function used in the construction of space-time planar point patterns. An stppp object is like a ppp object,
##' but with extra components for (1) a vector giving the time at whcih the event occurred and (2) a time-window over which observations occurred.
##' Observations are assumed to occur in the plane and the observation window is assumed not to change over time.
##'
##' @param P an object
##' @param ... additional arguments
##' @return method stppp 
##' @seealso \link{stppp}, \link{stppp.ppp}, \link{stppp.list}
##' @export 

stppp <- function(P,...){
    UseMethod("stppp")
}



##' stppp.ppp function
##'
##' Construct a space-time planar point pattern from a ppp object
##'
##' @method stppp ppp
##' @param P a spatstat ppp object
##' @param t a vector of length P$n
##' @param tlim a vector of length 2 specifying the observation time window
##' @param ... additional arguments
##' @return an object of class stppp
##' @seealso \link{stppp}, \link{stppp.list}
##' @export 

stppp.ppp <- function(P,t,tlim,...){
    if (length(t)!=P$n){
        stop("Vector t must have length P$n")
    }
    if (length(tlim)!=2){
        stop("tlim must be a vector of length exactly 2 and sorted in ascending order")
    }
    if (tlim[1] >= tlim[2]){
        stop("Vector tlim must be sorted in ascending order")
    }
    if (any(as.integer(t)>as.integer(tlim[2])) | any(as.integer(t)<as.integer(tlim[1]))){
        stop("There is a t such that as.integer(t) is outside range of as.integer(tlim). Please amend tlim so all observations fall in observation time window.")
    }
    P$t <- t
    P$tlim <- tlim
    class(P) <- c("stppp","ppp")
    return(P)
}



##' stppp.list function
##'
##' Construct a space-time planar point pattern from a list object
##'
##' @method stppp list
##' @param P list object containing $data, an (n x 3) matrix corresponding to (x,y,t) values; $tlim, a vector of length 2 givign the observation time window; and $window giving an owin spatial observation winow, see ?owin for more details 
##' @param ... additional arguments
##' @return an object of class stppp
##' @seealso \link{stppp}, \link{stppp.ppp}, 
##' @export

stppp.list <- function(P,...){
    tl <- P$tlim
        
    if(is.null(P$window)){
        stop("Require argument window")
    }
    if (is.null(P$tlim)){
        stop("Require argument tlim.")
    }

    idx <- inside.owin(P$data[,1],P$data[,2],P$window)
    if(sum(idx)!=length(idx)){
        warning(paste(sum(!idx),"points outside owin excluded"))
    }
    return(stppp(ppp(x=P$data[idx,1],y=P$data[idx,2],window=P$window),t=P$data[idx,3],tlim=P$tlim))            
    
}


##' print.stppp function
##'
##' Print method for stppp objects
##'
##' @method print stppp
##' @param x an object of class stppp
##' @param ... additional arguments
##' @return prints the stppp object x
##' @export 

print.stppp <- function(x,...){
    cat("Space-time point pattern\n")
    NextMethod("print",x)
    suppressWarnings(cat(paste("   Time Window : [",x$tlim[1],",",x$tlim[2],"]\n")))
}



##' plot.stppp function
##'
##' Plot method for stppp objects
##'
##' @method plot stppp
##' @param x an object of class stppp
##' @param ... additional arguments passed to plot
##' @return plots the stppp object x
##' @export 

plot.stppp <- function(x,...){
    plot(x$window,...)
    points(x$x,x$y,pch="+",col="red",cex=0.5)
}


##' Extract.stppp function
##'
##' extracting subsets of an stppp object.
##'
##' @name Extract.stppp
##' @aliases "[.stppp"
##' @param x an object of class stppp
##' @param subset the subset to extract 
##' @return extracts subset of an stppp object
##' @usage "x[subset]"
##' @examples
##' \dontrun{xyt <- lgcpSim()}
##' \dontrun{xyt}
##' \dontrun{xyt[xyt$t>0.5]}


"[.stppp" <- function(x,subset){
    x$x <- x$x[subset]
    x$y <- x$y[subset]
    x$t <- x$t[subset]
    x$n <- length(x$x)
    return(x)
}

##' as.ppp.stppp function
##'
##' Convert from stppp to ppp. Can be useful for data handling.
##'
##' @method as.ppp stppp
##' @param X an object of class stppp
##' @param ... additional arguments
##' @param fatal logical value, see details in generic ?as.ppp
##' @return a ppp object without observation times
##' @export 

as.ppp.stppp <- function(X,...,fatal=TRUE){
    return(ppp(x=X$x,y=X$y,window=X$window,marks=X$marks))
}

##' rescale.stppp function
##'
##' Rescale an stppp object. Similar to rescale.ppp
##'
##' @method rescale stppp
##' @param X an object of class stppp
##' @param s scale as in rescale.ppp: x and y coordinaes are scaled by 1/s
##' @param unitname parameter as defined in ?rescale
##' @return a ppp object without observation times
##' @export 

rescale.stppp <- function(X,s,unitname){
    t <- X$t
    tlim <- X$tlim
    xnew <- rescale.ppp(X,s,unitname)
    xnew <- stppp(xnew,tlim=tlim,t=t)
    return(xnew) 
}



##' integerise function
##'
##' Generic function for converting the time variable of an stppp object.
##'
##' @param obj an object
##' @param ... additional arguments
##' @return method integerise
##' @seealso \link{integerise.stppp}
##' @export

integerise <- function(obj,...){
    UseMethod("integerise")
}



##' integerise.stppp function
##'
##' Function for converting the times and time limits of an stppp object into integer values. Do this before estimating mu(t), and hence
##' before creating the temporalAtRisk object. Not taking this step is possible in lgcp, but can cause minor complications connected with the scaling of mu(t).
##'
##' @method integerise stppp
##' @param obj an stppp object
##' @param ... additional arguments
##' @return The stppp object, but with integerised times.
##' @export

integerise.stppp <- function(obj,...){
    attr(obj,"truetimes") <- obj$t # save previous values of time for possible use later
    if(all(as.integer(obj$t)==obj$t)){
        return(obj) # if times are already integer - valued, then assume everything is okay
    }
    else{
        obj$t <- as.integer(obj$t)
        if(as.integer(obj$tlim[2])==obj$tlim[2]){
            obj$tlim[2] <- obj$tlim[2] - 1
        }
        obj$tlim <- as.integer(obj$tlim)
        return(obj)        
    }
}
