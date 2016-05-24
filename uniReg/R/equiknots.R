equiknots <- function(a,b,g,k,coinc){
###########################################################################################
# Determines g+2k+2 knots for the spline basis. The inner knots lie equidistant in [a,b]. #
# If coinc=T, k knots are equal to each a and b, otherwise the outer knots are also equi- #
# distant beyond [a,b].                                                                   #
###########################################################################################
    if(!(is.vector(a,mode="numeric") && length(a)==1 && is.finite(a))){stop("a should be a finite numeric vector of length 1.")}
    if(!(is.vector(b,mode="numeric") && length(b)==1 && is.finite(b))){stop("b should be a finite numeric vector of length 1.")}
    if(!(g%%1==0 && g>=0 && is.finite(g))){stop("g should be a finite whole number >=0.")}
    if(!(k%%1==0 && k>=0 && is.finite(k))){stop("k should be a finite whole number >=0.")}
    if(!is.logical(coinc)){stop("coinc should be TRUE or FALSE")}
    
    inner <- seq(a,b,length.out=g+2)
    if(coinc){
        outera <- rep(a,k)
        outerb <- rep(b,k)
    }else{
        width <- (b-a)/(g+1)
        outera <- a - (k:1*width)
        outerb <- b + (1:k*width)
    }
    return(c(outera,inner,outerb))
}
