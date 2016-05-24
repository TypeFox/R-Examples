transM<-function(A,init,inter=100,stage.names =c("All grps",1:(ncol(A))),leg.room = 1.5,...){
m<-matrix(nrow=ncol(A),ncol=inter)
m[,1]<-init

for (i in 2:inter){
m[,i]<-A%*%m[,i-1]
}

n<-apply(m,2,sum)
plot(seq(1,inter),n,type="n",ylim=c(0,max(n)*leg.room),...)
lines(seq(1,inter),n,col=1,lwd=3)
for(i in 1:ncol(A))
lines(seq(1,inter),m[i,],col=i+1)
res<-matrix(nrow=3,ncol=inter,data=apply(m,2,function(x){x/sum (x)}),dimnames=list(stage.names[-1],seq(1,inter)))
legend("topright",legend=stage.names,lwd=c(3,rep(1,ncol(A))),col=seq(1,(ncol(A)+1)))
res
}