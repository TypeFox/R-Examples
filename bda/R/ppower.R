##  To compute the power for ANOVA-type test for proportions

ppower <- function(effect.size, ngroup, n0, n1, range.p, alpha=0.05, gridsize=512)
{
    stopifnot(gridsize>0)
    if(gridsize<10)
        warning("Gridsize is small")
    if(missing(range.p)){
        p0 <- seq(0,1,length=gridsize)
    }else{
        if(length(range.p)==1){
            p0 <- range.p;
            gridsize <- 1
        }else if(length(range.p)==2){
            stopifnot(range.p[1]>=0 && range.p[2] <=1)
            p0 <- seq(range.p[1], range.p[2], length=gridsize)
        }else{
            if(any(range.p < 0 | range.p > 1))
                stop("Invalid value(s) for 'p0'")
            p0 <- sort(range.p)
            gridsize <- length(p0)
        }
    }
    p1 <- p0 + effect.size
    p1size <- sum(p1>=0 | p1<=1)
    if(p1size<0)
        stop("'effect.size' invalid")

    n0 <- round(n0); n1 <- round(n1)
    stopifnot(n0 > 0)
    stopifnot(n1 > 0)
    stopifnot(alpha>0 && alpha <1)
    ngroup <- round(ngroup)
    stopifnot(ngroup > 1)
    k <- choose(ngroup, 2)
    alpha <- 1-(1-alpha)^(1/k) # reuse variable

    pwr <- .Fortran(.F_ppower,
                    as.double(p0), as.integer(gridsize),
                    as.double(effect.size), 
                    as.double(alpha),
                    as.integer(c(n0,n1)),
                    pwr=double(1))$pwr;
    pwr
}
