`plot.cca` <-
function(x,...){
   #Set the plotting parameters
   ncv<-length(x$corr)
   oldpar<-par(no.readonly=TRUE)
   on.exit(par(oldpar))
   par(mfrow=c(ceiling(sqrt(ncv)),ceiling(sqrt(ncv))))
   #Plot the data on each canonical variate
   for(i in 1:ncv){
      plot(x$canvarx[,i],x$canvary[,i],xlab="X",ylab="Y",main=paste("Canonical Variate Plot - Variate",i,sep=" "),...)
      abline(mean(x$canvary[,i],na.rm=TRUE)-x$corr[i]*mean(x$canvarx[,i],na.rm=TRUE),x$corr[i])
      text(mean(x$canvarx[,i],na.rm=TRUE),mean(x$canvary[,i],na.rm=TRUE),label=paste("r=",round(x$corr[i],digits=2),sep=""),pos=1,srt=180/pi*atan(x$corr[i]))    
   }
   #Redundancy plots
   par(mfrow=c(1,1),ask=TRUE)
   h<-cbind(x$xvrd,x$yvrd)        #Get the redundancy information
   h<-rbind(c(x$xrd,x$yrd),h)
   colnames(h)<-c("X Given Y","Y Given X")
   rownames(h)<-c("Total",paste("CV",1:ncv))
   barplot(h,beside=TRUE,main="Canonical Variate Redundancy Plot",ylim=c(0,1),ylab="Fraction of Variance Explained",legend.text=rownames(h),col=rainbow(ncv+1))
   #Plot the loadings on each canonical variate using the helio plot
   par(mfrow=c(ceiling(sqrt(ncv)),ceiling(sqrt(ncv))))
   for(i in 1:ncv){
      helio.plot(x,cv=i,main=paste("Structural Correlations for CV",i))
   }
   #Now repeat for variances deposited
   par(mfrow=c(ceiling(sqrt(ncv)),ceiling(sqrt(ncv))))
   for(i in 1:ncv){
      helio.plot(x,cv=i,main=paste("Explained Variance for CV",i),type="variance",axis.circ=c(0.5,1),range.rad=25)
   }
}

