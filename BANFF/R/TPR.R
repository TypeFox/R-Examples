####True Positive Rate
TPR=function(num,Trace)
{
  new1=Trace[,num]
  ratiosum=length(new1[which(new1>=1)])/(length(new1))
  
  finalratio=mean(ratiosum)
}
