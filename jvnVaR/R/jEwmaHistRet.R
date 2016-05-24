jEwmaHistRet <-
function(s, h){
# adjusted return by garch method
# initiating parameter
x0 <- 0.85
# estimating
 x <- jNewRapEwma(x0, s)
# volatility computing
nn <- length(s)
hh <- max(-nn+1,h)
vv <- s^2
s0 <- vv[1]
ss <- c()#vector of sigma^2 (volatility^2)
ss[1] <- s0
k <- min(nn+hh,nn+1)
for (i in 2:(nn+1)) {
ss[i] <- (1-x)*vv[i-1] + x*ss[i-1]
}
if (hh > 1){
for (i in k+1:nn+hh) {
ss[i] <- ss[i-1]
}
}
#adjus return
newRet <- s * sqrt(ss[nn+hh]) / sqrt(ss[1:nn])
return(newRet)
}
