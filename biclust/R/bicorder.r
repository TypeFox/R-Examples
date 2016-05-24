#### Function to order variables or objects that appear in a bicluster

bicorder<-function(bicResult, cols=TRUE, rev=FALSE)
{
le<-dim(bicResult@RowxNumber)[2]
res<-c()
if(!cols)
  {
  for(i in 1:le)
    {
    res<-c(res,which(bicResult@RowxNumber[,i])[!(which(bicResult@RowxNumber[,i]) %in% res)])
    }
  count<-1:dim(bicResult@RowxNumber)[1]
  res<-c(res,count[!(count %in% res)])
  }
else
  {
  for(i in 1:le)
    {
    res<-c(res,which(bicResult@NumberxCol[i,])[!(which(bicResult@NumberxCol[i,]) %in% res)])
    }
  count<-1:dim(bicResult@NumberxCol)[2]
  res<-c(res,count[!(count %in% res)])
  
  }

if(rev) res<-rev(res)

res
}

