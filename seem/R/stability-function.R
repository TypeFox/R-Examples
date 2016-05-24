# -------------------------
# stability
sta.plot <- function(A,x0,t){
 det <- det(A); Tr <- sum(diag(A))
 D <- Tr^2-4*det
 V <- eigen(A)$vectors
 L <- eigen(A)$values
 nt <-length(t) 
 X <- matrix(nrow=nt,ncol=2)

 p <- solve(V,x0)
 for(i in 1:nt)
  X[i,] <- p[1]*V[,1]*exp(L[1]*t[i])+ p[2]*V[,2]*exp(L[2]*t[i])
 matplot(t,X,type="l", col=1)
 legend("topright",legend=c("X1","X2"),lty=1:2,col=1,cex=0.7)
 mtext(paste("Det=",det,"Trace=",Tr,"Discr=", D),cex=0.7)
 mtext(paste("Eval1=",round(L[1],3),"Eval2=",round(L[2],3)),cex=0.6,line=-1)
}

phase.plot <- function(A,x0,t,long=20){
 det <- det(A); Tr <- sum(diag(A))
 D <- Tr^2-4*det
 V <- eigen(A)$vectors
 L <- eigen(A)$values
 X <- matrix(nrow=length(t),ncol=2)

 p <- solve(V,x0)
 for(i in 1:length(t))
 X[i,] <- p[1]*V[,1]*exp(L[1]*t[i])+ p[2]*V[,2]*exp(L[2]*t[i])
 plot(X[,1],X[,2],xlab="X1",ylab="X2",type="l")
 arrows(Re(X[10,1]),Re(X[10,2]),Re(X[long,1]),Re(X[long,2]),length=0.1,lwd=1.7)
 mtext(paste("Det=",det,"Trace=",Tr,"Discr=", D),cex=0.7)
 mtext(paste("Eval1=",round(L[1],3),"Eval2=",round(L[2],3)),cex=0.6,line=-1)
}


