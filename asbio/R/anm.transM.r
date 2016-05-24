anm.transM<-function(A,init,inter=100,stage.names =c("All grps",1:(ncol(A))),leg.room = 1.5,anim.interval=0.1,...){
A<-as.matrix(A)
m<-matrix(nrow=ncol(A),ncol=inter)
for (i in 2:inter){
m[,1]<-init
m[,i]<-A%*%m[,i-1]
}
n<-apply(m,2,sum)
for(i in 1:inter){
dev.hold()
fig<-plot(seq(1,inter),n,type="n",ylim=c(0,max(n)*leg.room),...);grid(fig)
lines(seq(1,i),n[1:i],col=1,lwd=3)
for(j in 1:ncol(A)){
points(seq(1:i),m[j,][1:i],col=1+j,type="l")
}
legend("topright",legend=stage.names,lwd=c(3,rep(1,ncol(A))),col=seq(1,(ncol(A)+1)),bg="white")
dev.flush()
Sys.sleep(anim.interval)
}
}

