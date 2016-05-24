## Julia Bischof
## 10-09-2015

plotClonesCopyNumber<-function(copyNumber=NULL, withOutliers=TRUE, color="gray", title=NULL, PDF=NULL,...){
  if(length(copyNumber)==0){
    stop("--> Copy number vector is missing")
  }else{
    copyNumber<-as.numeric(copyNumber)
  }
  if(withOutliers==FALSE){
    q<-quantile(copyNumber, probs = 0.75)
    copyNumber<-copyNumber[which(copyNumber<=q)]
  }
  
  if(length(PDF)>0){
    pdf(file = paste(PDF,"_Clone-copy-number.pdf",sep=""),width = 7,height = 7,pointsize = 18)
  }
  par(mar=c(6,5,4,2))
  boxplot(copyNumber,col=color,xlab=" ",ylab="Copy number",main=title, sub=if(withOutliers==FALSE){"(without outliers)"}else{""})
  if(length(PDF)>0){
    dev.off()
  }
}


