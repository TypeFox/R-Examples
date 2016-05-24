plot.BBSGoF <-
function(x, ...){

Adjusted.pvalues=x$Adjusted.pvalues
data=x$data

ss=seq(0.001,0.999,by=0.001)

plot(x$n.blocks,x$Tarone.pvalues,ylab="pvalue",ylim=c(0,1),xlim=c(x$kmin,x$kmax),xlab="blocks",main="Tarone's test");abline(h=0.05,col=2,lty=2)
rug(x$n.blocks)
par(ask=T) 
plot(x$n.blocks,x$cor,ylab="cor",xlim=c(x$kmin,x$kmax),xlab="blocks",main="Within-block correlation")
rug(x$n.blocks)
par(ask=T) 
plot(ss,dbeta(ss,x$beta.parameters[1],x$beta.parameters[2]),type="l",ylab="density",xlab="probability",main="Beta density"); abline(v=x$p,col=2,lty=2)
par(ask=T) 
plot(x$n.blocks,x$effects,ylab="effects",xlab="blocks",xlim=c(x$kmin,x$kmax),main="Decision plot",ylim=c(min(x$effects)-5,x$SGoF));abline(h=x$SGoF,col=2,lty=2)
rug(x$n.blocks)
if(x$adjusted.pvalues==TRUE){
par(ask=T) 
plot(sort(data),sort(Adjusted.pvalues),main="BB-SGoF Adjusted p-values",xlab="Unadjusted p-values",ylab="Adjusted p-values",xlim=c(0,1),ylim=c(0,1))

abline(0,1,lty=2,lwd=1.5)
par(ask=F)

}

}
