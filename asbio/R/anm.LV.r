#--------------------------- Lotka-Volterra Competition -------------------------#

anm.LVcomp<-function(n1,n2,r1,r2,K1,K2,a2.1,a1.2,time=seq(0,200),ylab="Abundance",xlab="Time",interval=0.1,...){
y<-xstart<-c(n1=n1,n2=n2)
pars<-c(r1=r1,r2=r2,K1=K1,K2=K2,a2.1=a2.1,a1.2=a1.2) 

pr<-as.list(pars)
func<-function(time=time,xstart=xstart,pars=pars){
    n1<-xstart[1]
    n2<-xstart[2]
    with(as.list(pars),{
    dn1<-r1*n1*((K1-n1-(a1.2*n2))/K1)
    dn2<-r2*n2*((K2-n2-(a2.1*n1))/K2)
    res<-list(c(dn1,dn2))
    })}

out<-as.data.frame(rk4(xstart,time,func,pars))
r1.lab<-bquote(paste(r[1],"=",.(pr$r1)));r2.lab<-bquote(paste(r[2],"=",.(pr$r2)))
K1.lab<-bquote(paste(K[1],"=",.(pr$K1)));K2.lab<-bquote(paste(K[2],"=",.(pr$K2)))
a21.lab<-bquote(paste(alpha[21],"=",.(pr$a2.1)));a12.lab<-bquote(paste(alpha[12],"=",.(pr$a1.2)))
layout(matrix(c(1,1,0,2,2,rep(3,20)),5,5,byrow=TRUE))
for(i in min(time):max(time)){
    dev.hold()
    par(mar=c(0,0,0,0))
    N1.lab<-bquote(paste(N[1],"=",.(round(out$n1[i+1],0))));N2.lab<-bquote(paste(N[2],"=",.(round(out$n2[i+1],0))));t.lab<-bquote(paste("Time = ",.(i)))
    plot(seq(0,1),type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
    legend("center", legend=c(as.expression(N1.lab),as.expression(N2.lab)),cex=1.75,title=as.expression(t.lab),ncol=2,bty="n")
    plot(seq(0,1),type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",cex=1.65)
    legend("center",legend=c(as.expression(r1.lab),as.expression(K1.lab),as.expression(a12.lab),as.expression(r2.lab),as.expression(K2.lab),
    as.expression(a21.lab)),cex=1.5,ncol=2,bty="n",title="Parameters")
    par(mar=c(5,4,2,2))
    plot(out$time,out$n1,type="n",xlab=xlab,ylab=ylab,ylim=c(min(out$n1,out$n2),max(out$n1,out$n2)),main="Lotka-Volterra Competition",...)
    grid()
    points(out$time[0:i+1],out$n1[0:i+1],type="l",lty=1,col=2)
    points(out$time[0:i+1],out$n2[0:i+1],type="l",lty=2,col=4)
    legend("topright",col=c(2,4),lty=c(1,2), legend=c("Sp. 1","Sp. 2"),bg="white")
 dev.flush()
 Sys.sleep(interval)
}
}





#----------------------Lotka-Volterra exploitation------------------------#

anm.LVexp<-function(nh,np,rh,con,p,d.p,time=seq(0,200),ylab="Abundance",xlab="Time",interval=0.1,circle=FALSE,...){
y=xstart=c(nh=nh,np=np);pars=c(rh=rh,con=con,p=p,d.p=d.p)

pr<-as.list(pars)
func<-function(time=time,x=xstart,pars=pars){
nh<-x[1]
np<-x[2]
with(as.list(pars),{
dnh<-rh*nh-(p*nh*np)
dnp<-con*p*nh*np-(d.p*np)
res<-list(c(dnh,dnp))
})}
out<-as.data.frame(rk4(xstart,time,func,pars))
rh.lab<-bquote(paste(r[h],"=",.(pr$rh)));c.lab<-bquote(paste(c,"=",.(pr$con)));p.lab<-bquote(paste(p,"=",.(pr$con)))
p.lab<-bquote(paste(p,"=",.(pr$p)));dp.lab<-bquote(paste(d[p],"=",.(pr$d.p)))
ccol<-rainbow(max(time)-min(time)+1)
layout(matrix(c(1,1,0,2,2,rep(3,20)),5,5,byrow=TRUE))
for(i in min(time):max(time)){
    dev.hold()
    par(mar=c(0,0,0,0))
    Nh.lab<-bquote(paste(Prey,"=",.(round(out$nh[i+1],0))));Np.lab<-bquote(paste(Pred,"=",.(round(out$np[i+1],0))));t.lab<-bquote(paste("Time = ",.(i)))
    plot(seq(0,1),type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
    legend("center", legend=c(as.expression(Nh.lab),as.expression(Np.lab)),cex=1.7,title=as.expression(t.lab),ncol=2,bty="n")
    plot(seq(0,1),type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",cex=1.65)
    legend("center",legend=c(as.expression(rh.lab),as.expression(p.lab),as.expression(c.lab),as.expression(dp.lab)),
    cex=1.5,ncol=2,bty="n",title="Parameters")
    par(mar=c(5,4,2,2))
    if(circle==FALSE){
    plot(out$time,out$nh,type="n",xlab=xlab,ylab=ylab,ylim=c(min(out$nh,out$np),max(out$nh,out$np)),main="Lotka-Volterra Exploitation",...)
    grid()
    points(out$time[0:i+1],out$nh[0:i+1],type="l",lty=1,col=2)
    points(out$time[0:i+1],out$np[0:i+1],type="l",lty=2,col=4)
    legend("topright",col=c(2,4),lty=c(1,2), legend=c("Prey","Pred"),bg="white")
    }
    if(circle==TRUE){
    plot(out$nh,out$np,type="n",xlab="Prey",ylab="Predator",main="Lotka-Volterra Exploitation",...)
    grid()
    points(out$nh[0:i+1],out$np[0:i+1],type="l",lwd=1.5,col=ccol[i+1])
    points(out$nh[i+1],out$np[i+1],pch=19)
        }
    dev.flush()
    Sys.sleep(interval)
}
}



