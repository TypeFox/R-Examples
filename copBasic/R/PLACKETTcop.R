"PLACKETTcop" <-
function(u, v, para=NULL, ...) {
    T <- para[1]
    if(is.null(para)) {
       warning("Empty para argument, need value on [0,Inf]")
       return()
    }
    if(T < 0) {
       warning("Theta < 0, invalid parameter")
       return()
    }
    if(T == 1)         return(u*v)    # the product copula
    if(T == 0)         return(W(u,v)) # lower copula bounds
    if(! is.finite(T)) return(M(u,v)) # upper copula bounds
    cop <- 1+(T-1)*(u+v)
    cop <- cop - sqrt(cop^2 - 4*u*v*T*(T-1))
    cop <- cop / (2*(T-1))
    return(cop)
}
