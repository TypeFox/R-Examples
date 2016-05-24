"twot.permutation" <-
function(x1=DAAG::two65$ambient, x2=DAAG::two65$heated, nsim=2000, plotit=TRUE){
# oldpar<-par(mar=par()$mar-c(1,0,1,0))
# on.exit(par(oldpar))
n1 <- length(x1)
n2<-length(x2)
n<-n1+n2
x<-c(x1,x2)
dbar <- mean(x2)-mean(x1)
z <- array(,nsim)
  for(i in 1:nsim){
    mn <- sample(n,n2,replace=FALSE)
    dbardash <- mean(x[mn]) - mean(x[-mn])
    z[i] <- dbardash
  }
pval <- (sum(z >= abs(dbar)) + sum (z <= -abs(dbar)))/nsim
if(plotit){plot(density(z), xlab="", main="", yaxs="i", ylim=c(0,0.08), cex.axis=0.8)
abline(v=dbar)
abline(v=-dbar, lty=2)
mtext(side=3,line=0.5, text=expression(bar(x[2])-bar(x[1])), at=dbar)
mtext(side=3,line=0.5, text=expression(-(bar(x[2])-bar(x[1]))), at=-dbar)}
print(signif(pval,3))
invisible()
}
