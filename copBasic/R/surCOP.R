"surCOP" <-
function(u, v, cop=NULL, para=NULL, exceedance=TRUE, ...) {
    #str(cop)
    if(! exceedance) {
      u <- 1 - u; v <- 1 - v
    }
    return(u + v - 1 + cop(1-u,1-v, para=para, ...))
}
