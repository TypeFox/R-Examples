AR.ForeB <-
function(x,b,h,e,p1)
{
n <- nrow(x)
p <- length(b)-1

b1 <- b[p+1]
b2 <- b[1:p]
rx <- x[(n-p+1):n,1]

for( i in 1:h)
{
    f <- b1 + sum( b2 * rx[length(rx):(length(rx)-p+1)]) + e[as.integer(runif(1, min=1, max=nrow(e)))]
    rx <- c(rx,f)
}
f <- rx[(p+1):(p+h)]
return(f)
}
