plot.TCA <-
function(x, pc=NULL, ...){   
   if(x$K==1)  cat("only one loading vector calculated")
   if(x$K>1){
      if(x$cov.input){
          cat("input is a covariance matrix, please input a raw data \n")
          plot(c(0:x$K),c(0,x$pev),type="b",xlab="PC number", ylab="Proportion of variance explained", main="PC number vs. Variance Explained")
      }    
      if(!x$cov.input){
          par(mfrow=c(1,2))
          plot(c(0:x$K),c(0,x$pev),type="b",xlab="PC number", ylab="Proportion of variance explained", main="PC number vs. Variance Explained")
          K=x$K
          if(is.null(pc)){
              plot(x$PC[,1],x$PC[,2],xlab=paste("PC",1,sep=""),ylab=paste("PC",2,sep=""), main="Principal Components", type="p", cex=0.5)                                   
          }else{
              plot(x$PC[,pc[1]],x$PC[,pc[2]],xlab=paste("PC",pc[1],sep=""),ylab=paste("PC",pc[2],sep=""), main="Principal Components",type="p", cex=0.5)
          }     
      }
  }    
}
