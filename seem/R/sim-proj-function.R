sim.proj <- function(t, x0, proj.mat, vars.plot, tunit, xunit, pdfout=F){
 nx= length(x0); nt =length(t)
 varlab=paste("St",as.character(c(1:nx)))
 
x <- matrix(nx*nt, ncol=nx, nrow=nt)
x.p <- x; xs <- x 
lambda <- t;xtot <- t

x[1,] <- x0
xtot[1] <- sum(x[1,])
lambda[1] <- 1
# proportions as percent
x.p[1,] <- 100*(x[1,]/xtot[1])
# scaled with  respect to last as percent
xs[1,] <- 100*(x[1,]/x[1,nx])

for(i in 2:nt) {
  x[i,]<- proj.mat%*%x[i-1,]
  xtot[i] <- sum(x[i,])
  lambda[i] <- xtot[i]/xtot[i-1]
  x.p[i,] <- 100*(x[i,]/xtot[i])
  xs[i,] <- 100*(x[i,]/x[i,nx])
} 

mat <- matrix(1:4,2,2,byrow=T)
nf <- layout(mat, widths=rep(7/2,2), heights=rep(7/2,2), TRUE)
par(mar=c(4,4,3,.5),xaxs="i",yaxs="i")

if(pdfout==T){
pdf("chp10/proj-mat-out.pdf")

matplot(t,x[,vars.plot], type="l", xlim=c(0,1.5*max(t)),lty=1:nx, col=1,
         xlab=paste("t",tunit), ylab=paste("X",xunit))
legend (1.01*max(t), max(x), legend=varlab[vars.plot], lty=1:nx,col=1)

matplot(t,x.p[,vars.plot], type="l", xlim=c(0,1.5*max(t)),lty=1:nx, col=1,
         xlab=paste("t",tunit), ylab="Prop (%)")
legend (1.01*max(t), max(x.p), legend=varlab[vars.plot], lty=1:nx,col=1)

plot(t,lambda, type="l", xlim=c(0,1.5*max(t)),lty=1:nx, col=1,
         xlab=paste("t",tunit), ylab="Lambda")

plot(t,xtot, type="l", xlim=c(0,1.5*max(t)),lty=1:nx, col=1,
         xlab=paste("t",tunit), ylab=paste("X Total",xunit))

dev.off()
} else{

matplot(t,x[,vars.plot], type="l", xlim=c(0,1.5*max(t)),lty=1:nx, col=1,
         xlab=paste("t",tunit), ylab=paste("X",xunit))
legend (1.01*max(t), max(x), legend=varlab[vars.plot], lty=1:nx,col=1)

matplot(t,x.p[,vars.plot], type="l", xlim=c(0,1.5*max(t)),lty=1:nx, col=1,
         xlab=paste("t",tunit), ylab="Prop (%)")
legend (1.01*max(t), max(x.p), legend=varlab[vars.plot], lty=1:nx,col=1)

plot(t,lambda, type="l", xlim=c(0,1.5*max(t)),lty=1:nx, col=1,
         xlab=paste("t",tunit), ylab="Lambda")

plot(t,xtot, type="l", xlim=c(0,1.5*max(t)),lty=1:nx, col=1,
         xlab=paste("t",tunit), ylab=paste("X Total",xunit))

}

x <- round(cbind(t,x),3)
xp <- round(cbind(t,x.p),3)
lambda <- round(cbind(t, lambda),3)
xs <- round(cbind(t,xs),3)
xtot <- round(cbind(t,xtot),3)

return(list(x=x, xp=xp, xs=xs,lambda=lambda, xtot=xtot))
}
