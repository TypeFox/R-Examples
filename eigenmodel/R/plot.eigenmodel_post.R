"plot.eigenmodel_post" <-
function(x,...) {

  n<-dim(x$Y)[1]
  R<-dim(x$L_postsamp)[2]
  p<-dim(x$b_postsamp)[2]

  np<- 1*(R>0) + 1*(p>0) 
  
  if(np>0) {

    if(np==1) { par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0)) }
    if(np==2) { par(mfrow=c(2,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0)) }

    
      if(R>0) {
 
      plot(c(1,dim(x$L_postsamp)[1]),range(x$L_postsamp),type="n",
           ylab=expression(lambda),xlab="sample") 
      abline(h=0,col="gray")

      L_ord<-t(apply(x$L_postsamp,1,sort))

      for(j in 1:R) {lines(L_ord[,j],col=j+1) }

      tmp<-eigen(x$ULU_postmean) 
      L<- diag(tmp$val[order(-abs(tmp$val))[seq(1,R,length=R)] ]/n,nrow=R)
      U<- tmp$vec[,order(-abs(tmp$val))[seq(1,R,length=R)],drop=FALSE ]*sqrt(n)
 
      if(R>1) { plot(U,type="n",xlab="",ylab="")  
                  abline(h=0,col="gray") ;  abline(v=0,col="gray") 
                  addlines(U,1*(x$Y>median(x$Y,na.rm=TRUE)))
                  #points(U,pch=16,cex=.8,col="blue")
                  text(U[,1],U[,2],labels=rownames(x$Y),cex=.8,col="blue")
               } 
 

              }

    if(p>0) {
    
      b_psamp<-x$b_postsamp 
      b_psamp<-t(t(b_psamp)*apply( apply(x$X,3,c),2,sd,na.rm=TRUE))
   
      plot(c(1,dim(b_psamp)[1]),range(b_psamp),type="n", 
           ylab=expression(beta),xlab="iteration")
      abline(h=0,col="gray")
      for(j in 1:p) {lines(b_psamp[,j],col=j+1) }

      plot(c(1-1/(p+1),p+1/(p+1)),range(b_psamp),type="n",ylab=expression(beta),
          xlab="",xaxt="n")
      axis(side=1,at=1:p,labels=dimnames(x$X)[[3]])
    
      abline(h=0,col="gray")
      qb<-apply(b_psamp,2,quantile,c(0.025,.5,.975))  
      for(j in 1:p) {
        segments(j,qb[1,j],j,qb[3,j],col=1+j) 
        points(j,qb[2,j]) 
                     } 
               }

             }
                             }

