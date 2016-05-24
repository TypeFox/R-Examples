signifz<-function(x, digits=6) {
    magnitude <- floor(log10(abs(x)))
    scale <- 10^(digits-magnitude-1)
    signs <- sign(x)
    ret <- x
    g0 <- which(signs>=0)
    ret[g0] <- floor(x[g0]*scale[g0])/scale[g0]
    l0 <- which(signs<0) 
    ret[l0] <- ceiling(x[l0]*scale[l0])/scale[l0]
    return(ret)
}
