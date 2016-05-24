fun.mean.convert <-
function (x, param, val = 0) 
{
    L2 <- x[2]
    L3 <- x[3]
    L4 <- x[4]
    r <- 1/L2 * ((L3 + 1)^-1 - (L4 + 1)^-1) + val
    if (param == "rs") {
        r <- -1/L2 * ((L3 + 1)^-1 - (L4 + 1)^-1) + val
    }
    rr <- c(r, L2, L3, L4)
    if (is.na(fun.theo.mv.gld(rr, param = param)[1])) {
        warning("No finite first moment, returning original input values")
        return(x)
    }
    return(rr)
}
