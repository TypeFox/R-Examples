GeoBoxPlot <- function(x, areas, ...) {
   ##
   ## Function to create the geographic boxplots of Willmott et al. (2007).
   ##
   ## Arguments:
   ##
   ## 'x' numeric vector of data to be box-plotted.
   ## 'areas' numeric vector of same length as 'x' giving the areas corresponding to each value of 'x'.
   ## '...' optional arguments to the 'boxplot' function.  The argument 'plot' is not allowed.
   ##
   ## Value: A plot is produced.  The values of the statistics going into the
   ##	boxplot's are returned invisibly.
   ##
   if( any( is.na(x)) | any(is.na( areas))) warning("GeoBoxPlot: missing values are not handled.")
   out <- boxplot( x, plot=FALSE, ...)
   ox <- order( x)
   a <- areas[ ox]
   a.frac <- a/sum( a)
   a2 <- cumsum( a.frac)
   x <- x[ox]
   n <- length( x)
   out$stats[,1] <- c( min(x), x[ min( (1:n)[ a2>0.25])], x[ min( (1:n)[ a2>0.5])], x[ min( (1:n)[ a2>0.75])], max( x))
   out <- bxp(out, ...)
   invisible( out)
}

KernelGradFUN <- function(x,ktype="LoG", nx=10, ny=12, sigma=1) return(kernel2dsmooth(x,kernel.type=ktype, nx=nx, ny=ny, sigma=sigma))

S1 <- function(x, ...) {
    UseMethod("S1", x)
} # end of 'S1' function.

S1.SpatialVx <- function(x, ..., xhat, gradFUN="KernelGradFUN", time.point=1, model=1) {
    a <- attributes(x)

    ## Begin: Get the data sets
    if(!missing(time.point) && !missing(model)) dat <- datagrabber(x, time.point=time.point, model=model)
    else if(!missing(time.point)) dat <- datagrabber(x, time.point=time.point)
    else if(!missing(model)) dat <- datagrabber(x, model=model)
    else dat <- datagrabber(x)

    X <- dat$X
    Xhat <- dat$Xhat
    ## End: Get the data sets

    res <- S1.default(x=X, ..., xhat=Xhat, gradFUN=gradFUN)

    if(length(a$data.name) == a$nforecast + 2) {
        dn <- a$data.name[-(1:2)]
        vxname <- a$data.name[2]
    } else {
        dn <- a$data.name[-1]
        vxname <- a$data.name[1]
    }
    if(!is.numeric(model)) model.num <- (1:a$nforecast)[dn == model]
    else model.num <- model

    names(res) <- paste(vxname, " vs ", dn[model.num], sep="")

    attr(res, "time.point") <- time.point
    attr(res, "model") <- model

    return(res)
} # end of 'S1.SpatialVx' function.

S1.default <- function(x, ..., xhat, gradFUN="KernelGradFUN") {
   Xgrad <- do.call(gradFUN, c(list(x),list(...)))
   Ygrad <- do.call(gradFUN, c(list(xhat),list(...)))
   denom <- sum(colSums(pmax(abs(Xgrad),abs(Ygrad),na.rm=TRUE),na.rm=TRUE),na.rm=TRUE)
   numer <- sum(colSums(abs(Ygrad - Xgrad),na.rm=TRUE),na.rm=TRUE)
   return(100*numer/denom)
} # end of 'S1.default' function.

ACC <- function(x, ...) {
    UseMethod("ACC", x)
} # end of 'ACC' function.

ACC.SpatialVx <- function(x, ..., xclim=NULL, xhatclim=NULL, time.point=1, model=1) {
    a <- attributes(x)

    ## Begin: Get the data sets
    if(!missing(time.point) && !missing(model)) dat <- datagrabber(x, time.point=time.point, model=model)
    else if(!missing(time.point)) dat <- datagrabber(x, time.point=time.point)
    else if(!missing(model)) dat <- datagrabber(x, model=model)
    else dat <- datagrabber(x)

    X <- dat$X
    Xhat <- dat$Xhat
    ## End: Get the data sets

    res <- ACC.default(x=X, ..., xhat=Xhat, xclim=xclim, xhatclim=xhatclim)

    if(length(a$data.name) == a$nforecast + 2) {
        dn <- a$data.name[-(1:2)]
        vxname <- a$data.name[2]
    } else {
        dn <- a$data.name[-1]
        vxname <- a$data.name[1]
    }
    if(!is.numeric(model)) model.num <- (1:a$nforecast)[dn == model]
    else model.num <- model

    names(res) <- paste(vxname, " vs ", dn[model.num], sep="")

    attr(res, "time.point") <- time.point
    attr(res, "model") <- model

    return(res)
} # end of 'ACC.SpatialVx' function.

ACC.default <- function(x, ..., xhat, xclim=NULL, xhatclim=NULL) {
   if(!is.null(xclim)) x <- x - xclim
   if(!is.null(xhatclim)) xhat <- xhat - xhatclim
   return(cor(c(x),c(xhat)))
} # end of 'ACC.default' function.
