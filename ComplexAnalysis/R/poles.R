poles <-
function(func,R,track.plot=TRUE,control=list(tol=10^-15,h=10^-12,ite=50,unit.concentric=2,sector=12,draw=1)){
control0<-list(tol=10^-15,h=10^-12,ite=50,unit.concentric=2,sector=12,draw=1)
control<-c(control0[setdiff(names(control0),names(control))],control)
if(track.plot){plot(R*exp(1i*seq(0,2*pi,0.001)),type="l",xlim=c(-R,R),ylim=c(-R,R),xlab="Real",ylab="Imaginary")}
f<-NULL
recip<-paste("f<-function(",names(formals(func)),"){0*1i+1/(",body2string(func),")}",collapse="",sep="")
eval(parse(text=recip))
zeros<-NULL;x0<-t(t(c(Inf,Inf)))
start<-NULL
if(R<1){R<-1}
R<-round(R)
r<-(R/control$unit.concentric)/sqrt(R/control$unit.concentric)
rI<-0
for(i in 1:(control$unit.concentric*R)){rI<-c(rI,r*sqrt(i))}
tI<-seq(0,2*pi,length.out=control$sector+1)
  for(i in 1:(length(rI)-1)){
  for(j in 1:(length(tI)-1)){
     if(track.plot){
         lines(c(0,Re(R*exp(1i*tI[j]))),c(0,Im(R*exp(1i*tI[j]))),lty=3,col="grey")
         z<-rI[i]*exp(1i*seq(0,2*pi,0.001))
         lines(Re(z),Im(z),lty=3,col="grey")
         }
     start<-c(start,runif(control$draw,rI[i],rI[i+1])*exp(runif(control$draw,tI[j],tI[j+1])*1i))
  }}
for(k in 1:length(start)){
 track<-NULL
 convergence<-FALSE
 x<-t(t(c(Re(start[k]),Im(start[k]))));
 J<-matrix(NA,nrow=2,ncol=2)
    for(i in 1:control$ite){
    track<-c(track,x[1]+x[2]*1i)
    Fx<-f(x[1]+x[2]*1i);Fx<-t(t(c(Re(Fx),Im(Fx))))
    #if(sqrt(sum(Fx^2))<control$tol) {convergence<-TRUE;break}
    if(sqrt(sum((x-x0)^2))<control$tol){convergence<-TRUE;break}
    if(x[1]>1 & x[2]>1){dx<-sqrt(control$h)*abs(x)}else{dx<-t(t(rep(sqrt(control$h),2)))}
    Fdx1<-f((x[1]+dx[1])+x[2]*1i)
    Fdx2<-f(x[1]+(x[2]+dx[2])*1i)
    Fdx1<-t(t(c(Re(Fdx1),Im(Fdx1))))
    Fdx2<-t(t(c(Re(Fdx2),Im(Fdx2))))
    J[,1]<-(Fdx1-Fx)/dx[1]
    J[,2]<-(Fdx2-Fx)/dx[2]
    tryCatch(delta<-solve(J)%*%((-1)*Fx), error=function(e) NULL)
    x0<-x
    x<-x+delta
    }
 if(convergence==TRUE){zeros<-c(zeros,x[1]+x[2]*1i);if(track.plot){lines(track,col=k%%6+2)}}
} #end of all start's
if(length(zeros)>0){result<-unique(round(zeros,abs(log10(control$tol))-1));points(result,pch=19)}else{result<-NULL}
if(track.plot){points(start,pch=21,cex=0.5)}
return(result)}
