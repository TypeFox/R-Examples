"onetPermutation" <-
function(x=DAAG::pair65$heated-DAAG::pair65$ambient, nsim=2000, plotit=TRUE){
oldpar<-par(mar=par()$mar-c(1,0,1,0))
on.exit(par(oldpar))
n <- length(x)
dbar <- mean(x)
absx <- abs(x)
z <- array(,nsim)
  for(i in 1:nsim){
    mn <- sample(c(-1,1),n,replace=TRUE)
    xbardash <- mean(mn*abs(x))
    z[i] <- xbardash
  }
pval <- (sum(z >= abs(dbar)) + sum (z <= -abs(dbar)))/nsim
if(plotit){plot(density(z), xlab="", main="", yaxs="i", cex.axis=0.8, bty="L")
abline(v=dbar)
abline(v=-dbar, lty=2)
mtext(side=3,line=0.5, text=expression(bar(d)), at=dbar)
mtext(side=3,line=0.5, text=expression(-bar(d)), at=-dbar)}
print(signif(pval,3))
invisible()
}
