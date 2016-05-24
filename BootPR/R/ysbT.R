ysbT <-
function(x,b,e)
{
    p <- length(b)-2
    n <- nrow(x)
    tm <- 1:n
    y <- x
    for(i in (n-p):1)
        y[i,1] <- b[p+1,1] + b[p+2,1]*tm[i] + sum(b[1:p,1]*y[(i+1):(i+p),1]) + e[i]
    return(y)
}
