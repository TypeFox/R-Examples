vimp.plot <-
function(fit, num=10, type=2, norm=TRUE, titles=TRUE)
{
 if (type%in%c(0,1,2)==FALSE) 
    stop("Specified plot type does not exist")
 tPIs<-sort(fit$PI.importance, decreasing=TRUE)[1:num]
 tpreds<-sort(fit$Predictor.importance, decreasing=TRUE)[1:num]
 if (norm==TRUE) {
     tPIs.norm<-tPIs/tPIs[1]
     tpreds.norm<-tpreds/tpreds[1]
     up1<-1.1
     up2<-1.1
     if (titles==TRUE) {titles<-c("Predictor Importance", "normalized scores", "PI importance")}
     else {titles<-c("","","")}  
     }
 else {
    tPIs.norm<-tPIs
    tpreds.norm<-tpreds
    up1<-tpreds[1]+.5
    up2<-tPIs[1]+.5
    if (titles==TRUE) {titles<-c("Predictor Importance", "variable importance scores", "PI importance")}
    else {titles<-c("","","")}
    }
 if (type==0)
  {
  par(cex.axis=0.6, mai=c(1,1.4,.5,.1), pch=19, las=1)
  plot(x<-tpreds.norm, y<-c(num:1), xlab=titles[2], ylab="", xlim=c(0.025, up1),
            ylim=c(0.5, num+.5), yaxt="n", main=titles[1])
  segments(0, num, tpreds.norm[1], num, lty=1)
  for (i in 1:(num-1))
    {
    segments(0,(num-i),tpreds.norm[i+1],(num-i),lty=1)
    } 
  axis(2, at=c(num:1), labels=names(tpreds))
  }
 if (type==1)
  {
  par(cex.axis=0.6, mai=c(1,1.4,.5,.1), pch=19, las=1)
  plot(x<-tPIs.norm, y<-c(num:1), xlab=titles[2], ylab="", xlim=c(0.025,up2),
            ylim=c(0.5, num+.5), yaxt="n", main=titles[3])
  segments(0, num, tPIs.norm[1], num, lty=1)
  for (i in 1:(num-1))
    {
    segments(0,(num-i),tPIs.norm[i+1],(num-i),lty=1)
    } 
  axis(2, at=c(num:1), labels=names(tPIs))
  }
 if (type==2)
  {
  par(mfrow=c(1,2), cex.axis=0.6, mai=c(1,1.4,.5,.1), pch=19, las=1)
  plot(x<-tpreds.norm, y<-c(num:1), xlab=titles[2], ylab="", xlim=c(0.025, up1),
            ylim=c(0.5, num+.5), yaxt="n", main=titles[1])
  segments(0, num, tpreds.norm[1], num, lty=1)
  for (i in 1:(num-1))
    {
    segments(0,(num-i),tpreds.norm[i+1],(num-i),lty=1)
    } 
  axis(2, at=c(num:1), labels=names(tpreds))
  plot(x<-tPIs.norm, y<-c(num:1), xlab=titles[2], ylab="", xlim=c(0.025,up2),
            ylim=c(0.5, num+.5), yaxt="n", main=titles[3])
  segments(0, num, tPIs.norm[1], num, lty=1)
  for (i in 1:(num-1))
    {
    segments(0,(num-i),tPIs.norm[i+1],(num-i),lty=1)
    } 
  axis(2, at=c(num:1), labels=names(tPIs))
  }
}
