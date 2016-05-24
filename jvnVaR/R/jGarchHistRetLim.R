jGarchHistRetLim <-
function(s,h,L,U){
# adjusted return by garch method
# initiating parameter 
x0 <- c(0.05,0.85,mean(s^2)*0.1)
# estimating
 x <- jNewRapGarchLim(x0, s,L,U)
# volatility computing
nn <- length(s)
hh <- max(-nn+1,h)
vv <- s^2
s0 <- x[3]/(1 - x[1] - x[2])
ss <- c() #vector of sigma^2 (volatility^2)
ss[1] <- s0
k <- min(nn+hh,nn+1)
for (i in 2:(nn+1)) {
ss[i] <- x[3] + x[1]*vv[i-1] + x[2]*ss[i-1]
}
if (hh > 1){
for (i in k+1:nn+hh) {
ss[i] <- x[3] + (x[1] + x[2])*ss[i-1]
}
}
#adjus return
newRet <- s * sqrt(ss[nn+hh]) / sqrt(ss[1:nn])
return(newRet)
}
