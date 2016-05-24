predplot=function(prior,n,yobs)
{
	y=0:n; a=prior[1]; b=prior[2]
	probs=pbetap(prior,n,y)
      m=max(probs)*1.05
	plot(y,probs,type="h",ylab="Probability",ylim=c(0,m),
 main=paste("Predictive Dist., beta(",a,",",b,") prior, n=",n,
            ", yobs=",yobs),lwd=2,col="blue")
	points(yobs,0,pch=19,cex=2.5,col="red")
      text(yobs,m/8,"yobs",col="red")}