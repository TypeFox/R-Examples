plotClass.result<-function(true.classified, cross.method, class.method, flag.feature, feat.num)
{

     if((cross.method=="sub-sampling")||(cross.method=="fold-crossval"))
    {
      if(flag.feature)
      {
        if(length(class.method)==1)
        {
          vrem=matrix(1:ncol(true.classified),ncol(true.classified),1)
        }
        else
        {
          vrem=sapply(class.method,function(z) which(colnames(true.classified)==z))
        }

        par(mfrow=c(length(class.method),1))
        for(i in 1:length(class.method))
        {
          if(feat.num==1)
          {
            values=true.classified[,vrem[i],drop=FALSE]
          }
          else
          {
            values=true.classified[,vrem[,i],drop=FALSE]
          }
          colnames(values)=seq(feat.num,1,-1)
          box=boxplot(values,main = paste("Cross validation",class.method[i]),ylab = "Classification accuracy",xlab="n of features",col=i+1)
        }
      }
      else
      {
        box=boxplot(true.classified,main = "Cross validation",ylab = "Classification accuracy",xlab="Classifiers",col=3)
      }
    }

    if(cross.method=="leaveOneOut")
    {
      if(flag.feature)
      {
        dim(true.classified)=c(length(class.method),feat.num)
        barplot (true.classified, beside=TRUE, col=(1:length(class.method))+1, border='white'
                 , xlab="n of features", ylab="Accuracy", names.arg=as.character(seq(feat.num,1,-1))
                 , ylim=c(0,1), main="Classification with n of features")
        legend("bottomright", col=(1:length(class.method))+1, class.method, bg="white", lwd=1, pch=1:length(class.method))

      }
      else
      {
        barplot (true.classified, col=(1:length(class.method))+1, border='white'
                 , space=0.2, xlab="Classifiers", ylab="Accuracy", names.arg=class.method
                 , ylim=c(0,1), main="Classification results")
        legend("bottomright", col=(1:length(class.method))+1, paste(class.method,format(true.classified,digits=3)),bg="white", lwd=1, pch=1:length(class.method))
      }
    }
}
