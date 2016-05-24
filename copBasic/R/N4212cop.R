"N4212cop" <-
function(u, v, para=NULL, infis=100, ...) {
    T <- para[1]

    if(is.null(para)) {
       warning("Empty para argument, need value on [1,Inf)")
       return()
    }

    if(T < 1) {
       warning("Theta < 1, invalid parameter")
       return()
    }

    if(T == 1)    return(PSP(u,v)) # the PSP copula
    if(T > infis) return(M(u,v)  ) # upper copula bounds

    cop <- ( u^-1 - 1 )^T + (v^-1 - 1)^T
    cop <- ( 1 + cop^(1/T) )^-1

    return(cop)
}

