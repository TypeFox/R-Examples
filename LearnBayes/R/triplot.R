triplot=function(prior,data,where="topright")
{
a=prior[1]; b=prior[2]
s=data[1]; f=data[2]

p = seq(0.005, 0.995, length = 500)
prior=dbeta(p,a,b)
like=dbeta(p,s+1,f+1)
post=dbeta(p,a+s, b+f)

m=max(c(prior,like,post))

plot(p,post,type="l", ylab="Density", lty=2, lwd=3,
 main=paste("Bayes Triplot, beta(",a,",",b,") prior, s=",s,", f=",f),
 ylim=c(0,m),col="red")
lines(p,like,lty=1, lwd=3,col="blue")
lines(p,prior,lty=3, lwd=3,col="green")
legend(where,c("Prior","Likelihood","Posterior"),
  lty=c(3,1,2), lwd=c(3,3,3), col=c("green","blue","red"))

}