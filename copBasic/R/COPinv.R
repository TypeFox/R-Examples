"COPinv" <-
function(cop=NULL, u, t, para=NULL, ...) {
    if(t == 0) return(0)
    if(t == 1) return(u)
    if(u == t) return(1)
    #str(cop)
    func <- function(x, u, LHS, cop, para=para, ...) {
            return(LHS - COP(cop=cop, u=u, v=x, para=para, ...))
    }
    lo <- COP(0,0, cop=cop, para=para, ...)
    if(is.na(lo) | ! is.finite(lo)) lo <- .Machine$double.eps
    my.rt <- NULL
    try(my.rt <- uniroot(func,interval=c(lo,1),
                         u=u, LHS=t, cop=cop, para=para, ...))

    if(is.null(my.rt)) return(NA)
    if(length(my.rt$root) != 0) {
      v <- my.rt$root
      return(v)
    } else {
      message(c("COPinv: ",u,t,"\n"))
      return(NA)
    }
}
