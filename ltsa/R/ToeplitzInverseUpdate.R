`ToeplitzInverseUpdate` <-
function( GI, r, rnew )
{   
    g <- rev(c(r,rnew)[-1])
    GIg <- c(crossprod(GI,g))
    e <- 1/(r[1]-sum(g*GIg))
    f <- -GIg*e
    A <- GI + outer(GIg,GIg)*e
    rbind(cbind(A,f),c(f,e))
}

