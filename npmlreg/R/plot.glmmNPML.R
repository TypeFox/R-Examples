"plot.glmmNPML" <-
function(x,plot.opt=15, noformat= FALSE, ...){

  k<-length(x$masses) 
  if (k==1){stop("No graphical output for k=1")} 

  if (!noformat){
      if  (plot.opt==0){
          stop("Specify plot.opt >0")
      } else if (plot.opt %in% c(1,2,4,8)){
          par(mfrow=c(1,1))
      } else if (plot.opt %in% c(3,5,6,9,10,12)){
          par(mfrow=c(2,1),cex=0.5,cex.axis=1.5,cex.lab=1.5)
      } else {
          par(mfrow=c(2,2),cex=0.5,cex.axis=1.1,cex.lab=1.1)
      }   
   }  
       
  if (plot.opt%%2==1){#Disparities
      plot(0:x$EMiter,x$Misc$Disparity.trend, col=1,type="l",xlab='EM iterations',ylab='-2logL')      
  }
  
  if (plot.opt %%4 %in% c(2,3) ){#  EM Trajectories
      R<- x$Misc$res; ylim<- x$Misc$ylim; followmass<-x$Misc$EMTraj     
      plot(0:x$EMiter,followmass[,1],col=1,type='l',ylim=ylim,ylab='mass points',xlab='EM iterations')
      for (i in 1:k){   
          lines(0:x$EMiter, followmass[,i],col=i)
          if (x$Misc$mform=="1"){points(rep(x$EMiter,length(R)),R)}
      }
  }
  
  if (plot.opt %%8 %in% c(4,5,6,7) ){#  EBP vs true values
       
      #class.col<- masspoint.classifier(x)
      class.col<- post(x, level="lower")$classif
      plot(x$y[1:length(x$weights)], predict(x,type="response"), xlab="true response", ylab="Emp. Bayes Pred." ,col=class.col)
      abline(0,1)         
  }
  
  if (plot.opt >7){   #wik
      if (is.infinite(abs(max(x$Misc$res)))){
          cat("Infinite values: Posterior probability plot not applicable for this object."); return()
      }
      pmax<-max(x$post.prob)
      plot(x$Misc$res, x$post.prob[,1],col=1,type="p",ylim=c(0,pmax),ylab="post.prob", xlab="Residuals")
      for (i in 2:k) {
          points(x$Misc$res, x$post.prob[,i],col=i,pch=i,type="p")
      }
  }
  
  invisible(x)
}
