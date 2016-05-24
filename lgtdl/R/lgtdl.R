lgtdl <- function(time, cov) {
  x<-data.frame(time=time,cov=cov)
  class(x) <- c("lgtdl", class(x))
 return(x)
 }

is.lgtdl <- function(x) if(inherits(x, "lgtdl")) TRUE else FALSE

as.lgtdl <- function(x, row.names = NULL) {
    if( is.lgtdl(x) ) return(x)
    mtrx<-is.matrix(x)
    if( !inherits(x, "data.frame") && !mtrx )
       stop("only data frames and matrices can be coerced to lgtdl")
    if( mtrx )
      x<-as.data.frame(x)
    dims <- dim(x)
    if(is.na(match("time", names(x))))
       names(x)[1]<- "time"
    class(x)<-c("lgtdl", class(x))
    if( any( diff(x$time) < 0 ) )
        x<-x[order(x$time), ]
    return(x)
}

getcov <- function(x, ...) UseMethod("getcov")

getcov.lgtdl <- function(x, cov, ...)
{
    nx<-names(x)
    tidx <- match("time", nx)
    if( is.na( tidx ) )
       stop("x must have a time component")
    if( missing(cov) || is.null(cov) ) {
       v<-1:length(nx)
       v<-v[-tidx]
       cov<-x[,v[1]]
    }
    else if( is.character(cov) || is.numeric(cov) )
       cov <- x[,cov]
    else
       stop("cov must be a character or a index")
    cov
}

covnames <- function(x)
{
    nx <- names(x)
    tidx <- match("time", nx)
    if( is.na(tidx) )
        stop("x must have a time component")
    return(nx[-tidx])
}

interpprev <- function(x, ...) UseMethod("interpprev")

interpprev.AsIs <- function(x, ...)
{
    lenx <- length(x)
    rval <- rep(NA, lenx)
    time <- list(...)[[1]]
    lt <- length(time)
    if( lt!= lenx ) {
        if(lt != 1)
            stop("time must be the same length as x")
        time <- rep(time,1)
    }
    for( i in 1:lenx )
        rval[i] <- interpprev(x[[i]], time[i])
    return(rval)
}


interpprev.lgtdl <- function(x, time, cov = NULL, ...)
{
    if( missing(time) )
     stop("time must be supplied")
    nx<-names(x)
    tidx <- match("time", nx)
    tvec <- x[,tidx]
    cov <- getcov(x, cov)
    ans <- approx(tvec, cov, time, "const")
    return( ans$y )
}

interplinear <- function(x, ...) UseMethod("interplinear")

interplinear.lgtdl <- function(x, time, cov = NULL, ...)
{
    if( missing(time) )
      stop("time must be supplied")
    nx <- names(x)
    tidx <- match("time", nx)
    tvec <- x[,tidx]
    cov <- getcov(x, cov)
    ans <- approx( tvec, cov, time, "linear")
    return( ans$y )
}

#when a column is extracted from a data frame its class is AsIs
#hence we need an AsIs method
interplinear.AsIs <- function(x,...)
{
        lenx <- length(x)
        rval <- rep(NA, lenx)
        time <- list(...)[[1]]
	lt <- length(time)
	if( lt!= lenx ) {
		if(lt != 1)
			stop("time must be the same length as x")
                time <- rep(time, lenx)
        }
        for( i in 1:lenx )
          rval[i] <- interplinear(x[[i]], time[i])
        return(rval)
}

plot.lgtdl <- function(x, ...)
{
    args<-list(...)
    xlab<-args$xlab
    if( is.null(xlab) )
        args$xlab <- "Time"
    cov <- getcov(x)
    dc <- dim(cov)
    if( is.null(dc) )
        ncov <- 1
    else
        ncov<- dc[2]
    cnames<- covnames(x)
    if( ncov > 1 ) {
        nrow <- ceiling(ncov/2)
        opar <- par(mfrow=c(2, nrow))
    }
    yl <- args$ylab
    args$x <- x[,"time"]
    if( is.null(args$type) )
	args$type <- "l"
    for( i in 1:ncov) {
        if( is.null(yl) )
           args$ylab<-cnames[i]
        else
            args$ylab <- yl[i]
        if( ncov == 1 )
            args$y <- cov
        else
            args$y <- cov[,i]
        do.call("plot", args)
    }
    if( ncov > 1 )
        par(opar)
}

#toString.lgtdl returns a length 1 character vector
#width is ignored
toString.lgtdl <- function(x, width, ...)
{
    n <- dim(x)[1]
    return(paste("lgtdl, length = ",n,sep=""))
}

print.lgtdl <- function(x, ...)
{
    oclass <- class(x)
    class(x) <- "data.frame"
    print(x)
    class(x) <- oclass
    invisible(x)
}
