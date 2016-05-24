ple.plot <-
function(dd,gam=0.95,p=0.95,xlow=0,xh=NA,ylow=0,yh=1, ...){
#
dname<- deparse(substitute(dd))

d<- plekm(dd,gam)
ylm<-c(ylow,yh)
if(is.na(xh)) xm<-max(d$a) else xm<-xh
xlm<-c(xlow,xm)
xcl<- unlist(percentile.ple(dd,p,gam,TRUE))
#
cv<- round(100*(1 - (1 -gam)*2))

t1<- paste( "PLE for ", dname," With  ",cv," Percent CLS \n",
    100*p,"th Percentile=   ",round(xcl[1],3),
    "   (",round(xcl[2],3),",",round(xcl[3],3),")",sep="" )


yl<-"PLE of F(x)"
tm<- stepfun( d$a,c(0,d$ple)) 
plot(tm ,xlim=xlm,ylim=ylm,pch=" ",col.points="red",main=t1,ylab=yl,...)
   points(d$a,d$ple,pch=21,bg="yellow")
 nx<-length(d$a);   x<- c(d$a,d$a[nx]) 
 y<-c(0,d$lower); y<- c(y[1:nx],y[nx],NA)
 lines(stepfun( x ,y ),col.hor="green4",lty=2,pch=" " )
 points(x[1:nx-1],y[2:nx],pch=25 ,bg="green2",ce=0.7)

x <- c(x[1],d$a) ; y<- d$ucbF;  y<- c(y[1],y);  y<- c(y,y[nx+1])
x <- d$a ; y<- c(NA,d$upper) ; y[length(y)]<-1
lines(stepfun( x,y ),col.hor="green4" ,lty=2,pch=" " )
      points(x[1:nx-1],y[2:nx],pch=24,bg="green2",cex=0.7)
abline(p,0,lty=1,col="red")
#abline(0.25,0,lty=9,col="gray25");abline(0.5,0,lty=9,col="gray25")

mle<- lnorm.ml(dd)
stp<- (xm - min(d$a)/10)/50
x<- seq(min(d$a)/10,xm,stp)
y<-plnorm(x,mle$mu,mle$sigma)
lines(x,y,col="blue",cex=0.5,lty=3,lwd=1.5)
invisible(d)
}

