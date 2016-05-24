
#---------------------------- Geometric Growth --------------------------#

anm.geo.growth<-function(n0,lambda,time=seq(0,20),ylab="Abundance",xlab="Time",interval=0.1,...){
n<-n0*lambda^time

layout(matrix(c(1,1,0,2,2,rep(3,20)),5,5,byrow=TRUE))
lambda.lab<-bquote(paste(lambda,"=",.(lambda)))
old.par <- par(no.readonly = TRUE)
for(i in min(time):max(time)){
    dev.hold()
    par(mar=c(0,0,0,0))
    N.lab<-bquote(paste(N,"=",.(round(n[i+1],0))));t.lab<-bquote(paste("Time = ",.(i)))
    plot(seq(0,1),type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
    legend("center", legend=c(as.expression(N.lab),as.expression(t.lab)),cex=1.75,bty="n")
    plot(seq(0,1),type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",cex=1.65)
    legend("center",legend=as.expression(lambda.lab),cex=1.5,ncol=2,bty="n")
    par(mar=c(5,4,2,2))
    plot(time,n,type="n",xlab=xlab,ylab=ylab,main="Geometric Population Growth",...)
    grid()
    points(time[0:i+1],n[0:i+1],type="l",lty=1,col=2)
 dev.flush()
 Sys.sleep(interval)
}
on.exit(par(old.par))
}

#--------------------------- Exponential Growth -------------------------#

anm.exp.growth<-function(n,rmax,time=seq(0,20),ylab="Abundance",xlab="Time",interval=0.1,...){
y<-xstart<-c(n=n)
pars<-c(rmax=rmax) 

pr<-as.list(pars)
func<-function(time=time,xstart=xstart,pars=pars){
    n<-xstart
    with(as.list(pars),{
    dn1<-rmax*n
    res<-list(dn1)
    })}

out<-as.data.frame(rk4(xstart,time,func,pars))
rmax.lab<-bquote(paste(r[max],"=",.(pr$rmax)))

layout(matrix(c(1,1,0,2,2,rep(3,20)),5,5,byrow=TRUE))
old.par <- par(no.readonly = TRUE)
for(i in min(time):max(time)){
    dev.hold()
    par(mar=c(0,0,0,0))
    N.lab<-bquote(paste(N,"=",.(round(out$n[i+1],0))));t.lab<-bquote(paste("Time = ",.(i)))
    plot(seq(0,1),type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
    legend("center", legend=c(as.expression(N.lab),as.expression(t.lab)),cex=1.75,bty="n")
    plot(seq(0,1),type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",cex=1.65)
    legend("center",legend=as.expression(rmax.lab),cex=1.5,ncol=2,bty="n")
    par(mar=c(5,4,2,2))
    plot(out$time,out$n,type="n",xlab=xlab,ylab=ylab,main="Exponential Population Growth",...)
    grid()
    points(out$time[0:i+1],out$n[0:i+1],type="l",lty=1,col=2)
 dev.flush()
 Sys.sleep(interval)
}
on.exit(par(old.par))
}

#--------------------------- Logistic Growth -------------------------#

anm.log.growth<-function(n,rmax,K,time=seq(0,60),ylab="Abundance",xlab="Time",interval=0.1,...){
y<-xstart<-c(n=n)
pars<-c(rmax=rmax,K=K) 

pr<-as.list(pars)
func<-function(time=time,xstart=xstart,pars=pars){
    n<-xstart
    with(as.list(pars),{
    dn1<-rmax*n*(1-(n/K))
    res<-list(dn1)
    })}

out<-as.data.frame(rk4(xstart,time,func,pars))
rmax.lab<-bquote(paste(r[max],"=",.(pr$rmax)))
K.lab<-bquote(paste(K,"=",.(pr$K)))
old.par <- par(no.readonly = TRUE)
layout(matrix(c(1,1,0,2,2,rep(3,20)),5,5,byrow=TRUE))
for(i in min(time):max(time)){
    dev.hold()
    par(mar=c(0,0,0,0))
    N.lab<-bquote(paste(N,"=",.(round(out$n[i+1],0))));t.lab<-bquote(paste("Time = ",.(i)))
    plot(seq(0,1),type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
    legend("center", legend=c(as.expression(N.lab),as.expression(t.lab)),cex=1.75,bty="n")
    plot(seq(0,1),type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",cex=1.65)
    legend("center",legend=c(as.expression(rmax.lab),as.expression(K.lab)),cex=1.5,ncol=2,bty="n",title="Parameters")
    par(mar=c(5,4,2,2))
    plot(out$time,out$n,type="n",xlab=xlab,ylab=ylab,main="Logistic Population Growth",...)
    grid()
    points(out$time[0:i+1],out$n[0:i+1],type="l",lty=1,col=2)
 dev.flush()
 Sys.sleep(interval)
}
on.exit(par(old.par))
}

