jGarch11Vol <-
function(s,x,h){
# Garch(1,1)
nn <- length(s)
hh <- max(-nn+1,h)
vv <- s^2
s0 <- x[3]/(1 - x[1] - x[2])
ss <- c()
# vector of sigma^2 (volatility^2)
ss[1] <- s0
if (nn+hh > 1){
k <- min(nn+hh,nn+1)
for (i in 2:k) {
ss[i] <- x[3] + x[1]*vv[i-1] + x[2]*ss[i-1]
}
if (hh > 1){
for (i in k+1:nn+hh) {
ss[i] <- x[3] + (x[1] + x[2])*ss[i-1]
}
}
}
object <- ss[nn+hh]
return(object)
}
