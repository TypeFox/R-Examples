"plot.ada" <-
function(x,kappa=FALSE,test=FALSE,cols= rainbow(dim(x$model$errs)[2]+1),tflag=TRUE,...){
  if(!inherits(x,"ada")){
    stop("Error:  Object is not of class ada")
  }
  kapstat<-function(tab=diag(2) ){
    if(dim(tab)[1]==1){
      return(0)
    }
    if(dim(tab)[1]==dim(tab)[2]){
      rs<-apply(tab,2,sum)
      cs<-apply(tab,1,sum)
      N<-sum(rs)
      E<-sum(rs*cs)/N^2
      O<-sum(diag(tab))/N
      return( (O-E)/(1-E) )
    }else{
      return(NA)
    }
  }
  mat<-x$model$errs
  iter<-x$iter
  k<-dim(mat)[2]/2
  if(length(cols)<(k*2))
    cols=rep(cols,k*2)
  old.par <- par(no.readonly = TRUE)
  if(kappa)
    op <- par(mfrow = c(1, 2),...)     
  odds<-1:k*2-1
  vals<-range(mat[,odds])
  vals[2]<-vals[2]*1.1
  plot(mat[,1],xlab=paste("Iteration",1,"to",iter),ylab="Error",
       ylim=vals,cex.lab=1,cex.main=1.2,type="l",col=cols[1],lwd=1)
  axis(1,at=iter,font=2)
  indx<-1:5 * floor((iter-5)/5)
  if(iter<=5){
    indx<-3
  }else{
    if(iter<=10){
      indx<-c(3,7)
    }else{
      if(iter<=15)
        indx<-c(3,7,12)
    }
  }
  points(indx,mat[indx,1],pch="1",cex=1.5)
  leg<-"Train"
  if(test & k>1){
     if(tflag)
       title("Training And Testing Error")
     matlines(1:iter,mat[,odds[-1]],col=cols[odds[-1]],lty=1)
     points(rep(indx,length(odds[-1])),as.vector(mat[indx,odds[-1]]),
            pch=paste(sort(rep(2:k,length(indx)))),cex=1.5)
     leg=c("Train",paste("Test",1:(k-1),sep=""))
  }else{
    if(tflag)
      title("Training Error")
  }
  legend(par()$usr[1],par()$usr[4],leg,pch=paste(1:k))  
  if(kappa){
    odds<-odds+1
    mat[,odds]=1-mat[,odds]
    vals<-range(mat[,odds])
    
    plot(mat[,2],xlab=paste("Iteration",1,"to",iter),ylab="Kappa Accuracy",
         ylim=vals,cex.lab=1,cex.main=1.2,type="l",col=cols[2],lwd=1)
    axis(1,at=iter,font=2)
    points(indx,mat[indx,2],pch="1",cex=1.5)
    leg="Train"
    if(test & k>1){
      if(tflag)
        title("Training And Testing Kappas")
      matlines(1:iter,mat[,odds[-1]],col=cols[odds[-1]],lty=1)
      points(rep(indx,length(odds[-1])),as.vector(mat[indx,odds[-1]]),
             pch=paste(sort(rep(2:k,length(indx)))),cex=1.5)
      leg=c("Train",paste("Test",1:(k-1),sep=""))
    }else{
      if(tflag)
        title("Training Kappa")
    }
      legend(par()$usr[1],par()$usr[4],leg,pch=paste(1:k))  
  }
  on.exit(par(old.par))
}
