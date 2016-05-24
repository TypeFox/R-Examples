SimulateHReg <-
function(n, f, A, B, simPlot=FALSE, Dist="n", phi=0.3){
t<-1:n
x<-A*cos(2*pi*f*t)+B*sin(2*pi*f*t)
e<-switch(Dist,n=rnorm(n), t=rt(n,5)/1.29099, s=(rt(n,5,6)-7.1365)/3.27675, a=SimulateAR1(n,phi))
z<-x+e
amp<-A^2 + B^2
if (simPlot) {
    ymax <- max(z, amp)
    ymin <- min(z, -amp)
    plot(t,z, xlab="t", ylab="z[t]", pch=16, cex=1.5, lwd=3,ylim=c(ymax,ymin))
    points(t,x,pch=16, col="blue",cex=1)
    tt<-seq(1,n,length.out=100)
    xx<-A*cos(2*pi*f*tt)+B*sin(2*pi*f*tt) 
    lines(tt,xx,col="blue")
    }
z
}

