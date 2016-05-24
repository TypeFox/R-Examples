"COPinv2" <-
function(cop=NULL, v, t, para=NULL, ...) {
    if(t == 0) return(0)
    if(t == 1) return(v)
    if(v == t) return(1)
    #str(cop)
    func <- function(x, v, LHS, cop, para=para, ...) {
            return(LHS - COP(cop=cop, v=v, u=x, para=para, ...))
    }
    lo <- COP(0,0, cop=cop, para=para, ...)
    if(is.na(lo) | ! is.finite(lo)) lo <- .Machine$double.eps
    my.rt <- NULL
    try(my.rt <- uniroot(func,interval=c(lo,1),
                          v=v, LHS=t, cop=cop, para=para, ...))

    if(is.null(my.rt)) return(NA)
    if(length(my.rt$root) != 0) {
      u <- my.rt$root
      return(u)
    } else {
      message(c("COPinv2: ",v,t,"\n"))
      return(NA)
    }
}
