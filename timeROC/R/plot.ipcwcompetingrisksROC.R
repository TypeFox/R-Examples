plot.ipcwcompetingrisksROC <- function(x,FP=2,time,col="red",add=FALSE,title=TRUE,...){  
  if (FP!=1 & FP!=2){
    stop("FP indicates the type of specificity. It must be 1 or 2 \n") }
  if (missing(time)){
    stop("time is missing. \n") }
  if (length(time)!=1){
    stop("time must be an unique value. \n It indicates the time point at which the ROC curve is plotted. \n") }
  if ( ! time %in% x$times){
    stop("time must be one of the value included in the vector times that you used  for estimation.\n") }  
  if (add==FALSE) {
    plot(0,0,type="l",xlim=c(0,1),ylim=c(0,1),xlab="1-Specificity",ylab="Sensitivity")
    if(title==TRUE){
      if(FP==2){title(paste("ROC at time t=",time,", AUC=",round(x$AUC_2[which(x$times==time)]*100,1),sep=""))}
      if(FP==1){title(paste("ROC at time t=",time,", AUC=",round(x$AUC_1[which(x$times==time)]*100,1),sep=""))}
    }
  }
  if(FP==2){
    lines(x$FP_2[,which(x$times==time)],x$TP[,which(x$times==time)],col=col,type="l",...)
  }
  if(FP==1){
    lines(x$FP_1[,which(x$times==time)],x$TP[,which(x$times==time)],col=col,type="l",...)
  }
  abline(0,1,lty=2)
}
