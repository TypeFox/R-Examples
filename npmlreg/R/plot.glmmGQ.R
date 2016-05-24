"plot.glmmGQ" <-
function(x,plot.opt=3, noformat=FALSE,  ...){

  k<-length(x$masses)  

  if(!noformat){
      if (k==1){
          stop("No graphical output for k=1")
      }
      if  (plot.opt==0){
          stop("Specify plot.opt >0")
      } else if (plot.opt %in% c(1,2)){
          par(mfrow=c(1,1))
      } else if (plot.opt ==3){
          par(mfrow=c(2,1),cex=0.5,cex.axis=1.5,cex.lab=1.5)
      }
  }
      
  if (plot.opt%%2==1){#Disparities
      plot(0:x$EMiter,x$Misc$Disparity.trend, col=1,type="l",xlab='EM iterations',ylab='-2logL')      
  }
  
  if (plot.opt %%4 %in% c(2,3) ){#  EBP vs true values #klappt
      class.col<-post(x, level="lower")$classif
      plot(x$y[1:length(x$weights)], predict(x,type="response"), xlab="true response", ylab="Emp. Bayes Pred." ,col=class.col)
      abline(0,1)         
  }
    
  invisible(x)
}

