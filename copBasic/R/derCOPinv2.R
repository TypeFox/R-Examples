"derCOPinv2" <-
function(cop=NULL, v, t,
         delv=.Machine$double.eps^0.50, para=NULL, ...) {

    func <- function(x,v,LHS,cop,delv=delv,para=para, ...) {
            LHS - derCOP2(cop=cop, u=x, v=v, delv=delv, para=para, ...)
    }
    my.rt <- NULL
    try(my.rt <- uniroot(func,interval=c(0,1), v=v, LHS=t,
                              cop=cop, delv=delv, para=para, ...))
    if(is.null(my.rt)) return(NA) # Now the returned root is "u"
    ifelse(length(my.rt$root) != 0, return(my.rt$root), return(NA))
}
