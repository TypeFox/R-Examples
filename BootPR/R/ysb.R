ysb <-
function(x,b,e)
{
    p <- length(b)-1
    n <- nrow(x)
    y <- x
    for(i in (n-p):1)
        y[i,1] <- b[p+1,1] + sum(b[1:p,1]*y[(i+1):(i+p),1]) + e[i]
    return(y)
}
