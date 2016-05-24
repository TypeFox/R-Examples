###
# mstppp class definition and functions
###

##' mstppp function
##'
##' Generic function used in the construction of marked space-time planar point patterns. An mstppp object is like an stppp object,
##' but with an extra component containing a data frame (the mark information).
##'
##' Observations are assumed to occur in the plane and the observation window is assumed not to change over time.
##'
##' @param P an object
##' @param ... additional arguments
##' @return method mstppp 
##' @seealso \link{mstppp}, \link{mstppp.ppp}, \link{mstppp.list}
##' @export 

mstppp <- function(P,...){
    UseMethod("mstppp")
}



##' mstppp.ppp function
##'
##' Construct a marked space-time planar point pattern from a ppp object
##'
##' @method mstppp ppp
##' @param P a spatstat ppp object
##' @param t a vector of length P$n
##' @param tlim a vector of length 2 specifying the observation time window
##' @param data a data frame containing the collection of marks
##' @param ... additional arguments
##' @return an object of class mstppp
##' @seealso \link{mstppp}, \link{mstppp.list}
##' @export 

mstppp.ppp <- function(P,t,tlim,data,...){
    
    if (length(t)!=P$n){
        stop("Vector t must have length P$n")
    }
    if(!is.data.frame(data)){
        stop("data must be an object of class data.frame")
    }
    if (nrow(data)!=P$n){
        stop("Data frame  must have P$n rows")
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
    P$data <- data
    class(P) <- c("mstppp","stppp","ppp")
    return(P)
}

##' mstppp.stppp function
##'
##' Construct a marked space-time planar point pattern from an stppp object
##'
##' @method mstppp stppp
##' @param P an lgcp stppp object
##' @param data a data frame containing the collection of marks
##' @param ... additional arguments
##' @return an object of class mstppp
##' @seealso \link{mstppp}, \link{mstppp.list}
##' @export 

mstppp.stppp <- function(P,data,...){  
    if(!is.data.frame(data)){
        stop("data must be an object of class data.frame")
    }  
    if (nrow(data)!=P$n){
        stop("Data frame  must have P$n rows")
    }
    P$data <- data
    class(P) <- c("mstppp","stppp","ppp")
    return(P)
}



##' mstppp.list function
##'
##' Construct a marked space-time planar point pattern from a list object
##'
##' @method mstppp list
##' @param P list object containing $xyt, an (n x 3) matrix corresponding to (x,y,t) values; $tlim, a vector of length 2 givign the observation time window, 
##' $window giving an owin spatial observation winow, see ?owin for more details, and $data, a data frame containing the collection of marks 
##' @param ... additional arguments
##' @return an object of class mstppp
##' @seealso \link{mstppp}, \link{mstppp.ppp}, 
##' @export

mstppp.list <- function(P,...){
    tl <- P$tlim
        
    if(is.null(P$window)){
        stop("Require argument window")
    }
    if (is.null(P$tlim)){
        stop("Require argument tlim.")
    }
    if (is.null(P$data)){
        stop("Require argument data.")
    }

    idx <- inside.owin(P$xyt[,1],P$xyt[,2],P$window)
    if(sum(idx)!=length(idx)){
        warning(paste(sum(!idx),"points outside owin excluded"))
    }
    return(mstppp(ppp(x=P$xyt[idx,1],y=P$xyt[idx,2],window=P$window),t=P$xyt[idx,3],tlim=P$tlim,data=P$data))            
    
}


##' print.mstppp function
##'
##' Print method for mstppp objects
##'
##' @method print mstppp
##' @param x an object of class mstppp
##' @param ... additional arguments
##' @return prints the mstppp object x
##' @export 

print.mstppp <- function(x,...){
    cat("Marked Space-time point pattern\n")
    NextMethod("print",x)
    cat(paste("Marks:\n $data: a data frame with",ncol(x$data),"columns","\n"))
}



##' plot.mstppp function
##'
##' Plot method for mstppp objects
##'
##' @method plot mstppp
##' @param x an object of class mstppp
##' @param cols optional vector of colours to plot points with
##' @param ... additional arguments passed to plot
##' @return plots the mstppp object x
##' @export 

plot.mstppp <- function(x,cols="red",...){
    plot(x$window,...)
    points(x$x,x$y,pch="+",col=cols,cex=0.5)
}


##' Extract.mstppp function
##'
##' extracting subsets of an mstppp object.
##'
##' @name Extract.mstppp
##' @aliases "[.mstppp"
##' @param x an object of class mstppp
##' @param subset subsetto extract
##' @usage "x[subset]"
##' @return extracts subset of an mstppp object


"[.mstppp" <- function(x,subset){
    x$x <- x$x[subset]
    x$y <- x$y[subset]
    x$t <- x$t[subset]
    x$data <- x$data[subset,]
    x$n <- length(x$x)
    return(x)
}

##' as.ppp.mstppp function
##'
##' Convert from mstppp to ppp. Can be useful for data handling.
##'
##' @method as.ppp mstppp
##' @param X an object of class mstppp
##' @param ... additional arguments
##' @param fatal logical value, see details in generic ?as.ppp
##' @return a ppp object without observation times
##' @export 

as.ppp.mstppp <- function(X,...,fatal=TRUE){
    return(ppp(x=X$x,y=X$y,window=X$window,marks=X$data))
}

##' rescale.mstppp function
##'
##' Rescale an mstppp object. Similar to rescale.ppp
##'
##' @method rescale mstppp
##' @param X an object of class mstppp
##' @param s scale as in rescale.ppp: x and y coordinaes are scaled by 1/s
##' @param unitname parameter as defined in ?rescale
##' @return a ppp object without observation times
##' @export 

rescale.mstppp <- function(X,s,unitname){
    t <- X$t
    tlim <- X$tlim
    data <- X$data
    xnew <- rescale.ppp(X,s,unitname)
    xnew <- mstppp(xnew,tlim=tlim,t=t,data=data)
    return(xnew) 
}



##' integerise.mstppp function
##'
##' Function for converting the times and time limits of an mstppp object into integer values.
##'
##' @method integerise mstppp
##' @param obj an mstppp object
##' @param ... additional arguments
##' @return The mstppp object, but with integerised times.
##' @export

integerise.mstppp <- function(obj,...){
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
