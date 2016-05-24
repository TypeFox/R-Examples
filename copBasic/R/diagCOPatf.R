"diagCOPinv" <- function(f, cop=NULL, para=NULL, verbose=FALSE, interval=NULL,
                                  tol=.Machine$double.eps/10, ...) {
   diagCOPatf(f,cop=cop,para=para,verbose=verbose,interval=interval,tol=tol,...)
}

"diagCOPatf" <-
 function(f, cop=NULL, para=NULL, verbose=FALSE, interval=NULL,
                                  tol=.Machine$double.eps/10, ...) {
    if(verbose & length(f) > 1) {
       warning("f argument is not a single value, only first value will be used, ",
               "because verbose return needed")
       f <- f[1]
    }
    if(is.null(cop)) {
       warning("cop argument is NULL, returning NULL")
       return(NULL)
    }

    intv <- c(.Machine$double.eps, 1-.Machine$double.eps)
    if(! is.null(interval)) intv <- interval
    if(intv[1] < 0) {
       warning("invalid first value for the interval, returning NULL")
       return(NULL)
    }
    if(intv[2] > 1) {
       warning("invalid second value for the interval, returning NULL")
       return(NULL)
    }

    "afunc" <- function(uv, aF=NA) (aF - cop(uv, uv, para=para, ...))

    DatF <- sapply(f, function(af) {
                   if(af == 0) return(0)
                   if(af == 1) return(1)
                   rt <- NULL
                   try(rt <- uniroot(afunc, interval=intv, tol=tol, aF=af, ...))
                   if(is.null(rt)) {
                      warning("could not root for the diagonal inverse")
                      return(NA)
                   }
                   if(verbose) return(rt)
                   return(rt$root)
            })
    return(DatF)
 }

