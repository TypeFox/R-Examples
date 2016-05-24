.check <- function(x, sign='none',inf.rm=TRUE){
    stopifnot(is.numeric(x))
    sign <- match.arg(tolower(sign),c('nn','np','none'))
    size <- length(x)
    x.na <- is.na(x)
    n.na <- sum(x.na)
    if(any(x.na)){
        x <- x[!x.na]
        warning(paste(n.na, "'NA' value(s) removed"))
    }
    x.finite <- is.finite(x)
    n.inf <- sum(!x.finite)
    if(any(!x.finite)&&inf.rm){
        x <- x[x.finite]
        warning(paste(n.inf, "'Inf' value(s) removed"))
    }
    if(sign=='nn'){
        if(any(x<0))
            stop("Negative value not allowed")
    }else if(sign=='np'){
        if(any(x > 0))
            stop("Positive value not allowed")
    }
    list(x=x, size=size, n.inf=n.inf, n.na=n.na)
}

.edf <- function(x){
    y <- .check(x, inf.rm=FALSE)$x
    y.t <- table(y)
    Fy <- cumsum(y.t)/sum(y.t)
    list(x=as.numeric(names(y.t)), y=Fy)
}

.onUnload <- function(libpath)
  library.dynam.unload("bda",  libpath)

.bdaConnect <- NULL

