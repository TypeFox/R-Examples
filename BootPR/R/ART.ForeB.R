ART.ForeB <-
function(x,b,h,e,p1)
{
n <- nrow(x)
p <- length(b)-2
tm <- (n+1):(n+h)

b1 <- b[p+1]
b2 <- b[p+2]
b3 <- b[1:p]
rx <- x[(n-p+1):n,1]

for( i in 1:h)
{
    f <- b1 + b2*tm[i] + sum( b3 * rx[length(rx):(length(rx)-p+1)]) +e[as.integer(runif(1, min=1, max=nrow(e)))]
    rx <- c(rx,f)
}
f <- rx[(p+1):(p+h)]
return(f)
}
