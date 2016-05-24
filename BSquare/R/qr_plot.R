qr_plot <-
function(object,index=1,
                  xlab="Quantile level",ylab="Covariate effect",main="",
                  col=gray(.75),lwd=1,add=FALSE){
  tau<-object$tau
  q<-apply(object$q[,,index],2,quantile,c(0.05,0.5,0.95))
  ylim<-range(q)+0.25*(max(q)-min(q))*c(-1,1)
  
  if(!add){
    plot(NA,xlim=range(tau),ylim=ylim,xlab=xlab,ylab=ylab,main=main)
  }
  
  L<-q[1,]
  U<-q[3,]
  for(l in 1:(length(tau)-1)){      
    xxx<-tau[l+c(0,0,1,1,0)] 
    yyy<-c(L[l],U[l],U[l+1],L[l+1],L[l])
    polygon(xxx,yyy,col=col,border=FALSE)
  }
  lines(tau,q[1,])
  lines(tau,q[2,],lwd=1.5*lwd)
  lines(tau,q[3,])
}
