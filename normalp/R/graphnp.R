graphnp<-function(p=c(1,2,3,4,5),mu=0,sigmap=1,title="Exponential Power Distributions"){
if(!is.numeric(p)||!is.numeric(mu)||!is.numeric(sigmap))
stop (" Non-numeric argument to mathematical function")
n<-length(p)
leg<-vector()
for (i in 1:n){
if (p[i]<1) stop("There is a value of p less than 1")
}
if (n>5) stop("Please, insert only five values for p")
n1<-0
xfit<-seq(mu-4*sigmap,mu+4*sigmap,0.01)
n1<-length(p[p>=50])
if (n==1 && p>=50 || n==n1){
plot(xfit,dunif(xfit,min=mu-sigmap,max=mu+sigmap),xlim=range(mu-4*sigmap,mu+4*sigmap),
ylim=range(0,(1+0.1)/(2*sigmap)),
type="l",main=title,col=1,xlab="x",ylab="f(x)")
leg<-c(leg,expression(paste("p-> ",infinity)))
legend(mu+2.5*sigmap, (1-0.1)/(2*sigmap), leg,col=1, lty=1)
}
else {
p<-sort(p)
n<-n-n1
p<-p[1:n]
if (n1>0) { plot(xfit,dunif(xfit,min=mu-sigmap,max=mu+sigmap),xlim=range(mu-4*sigmap,mu+4*sigmap),ylim=range(0,(1+0.1)/(2*sigmap)),type="l",main=title,col=1,xlab="x",ylab="f(x)")
        if(n>0) for(i in 1:n){lines(spline(xfit,dnormp(xfit,mu=mu,sigmap=sigmap,p=p[i])),col=i+1)}
          leg<-paste("p=",p)
        leg<-c(leg,expression(paste("p-> ",infinity)))
        legend(mu+2.5*sigmap, (1-0.1)/(2*sigmap), leg,col=c(1:n+1,1), lty=1)
}
else {
    plot(xfit,dnormp(xfit,mu=mu,sigmap=sigmap,p=p[1]),xlim=range(mu-4*sigmap,mu+4*sigmap),ylim=range(0,(1+0.1)/(2*sigmap)),type="l",main=title,col=1,xlab="x",ylab="f(x)")
        if(n>1) for(i in 2:n){lines(spline(xfit,dnormp(xfit,mu=mu,sigmap=sigmap,p=p[i])),col=i)}
    leg<-paste("p=",p)
    legend(mu+2.5*sigmap, (1-0.1)/(2*sigmap), leg,col=1:n, lty=1)
    }
} 
}

