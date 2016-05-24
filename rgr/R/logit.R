logit <-
function (pp) 
{
    temp.p <- remove.na(pp)
    p <- temp.p$x[1:temp.p$n]
    for (i in 1:temp.p$n) {
        if ((p[i] < 0) | (p[i] > 1)) 
            stop("The proportion must be in the range 0 to 1")
    }
    z <- log(p) - log(1 - p)
    return(z = z)
}
