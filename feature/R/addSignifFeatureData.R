########## R function: addSignifFeatureData ##########

# For adding data points which lie in significant features
# regions to a KDE plot.

# Last changed: 04 NOV 2005

addSignifFeatureData <- function(x,SignifFeatureData,featureCol,dest,
   trans.alpha=1)
{
  if (is.vector(x))
    d <- 1
  else
    d <- ncol(x)
  
  if (d==1)
    if (sum(SignifFeatureData)>0)
      rug(x[SignifFeatureData], col=featureCol)
      
  if (d==2)
    points(x[SignifFeatureData,],pch=1,col=featureCol) 
  
  if (d==3)
  {
    if (!all(SignifFeatureData==FALSE))
    {
      x.sig <- x[SignifFeatureData,]
      if (is.vector(x.sig))
        x.sig <- matrix(x.sig, nrow=1)
      
      points3d(x.sig[,1],x.sig[,2], x.sig[,3], size=3,color=featureCol, alpha=trans.alpha)           
    }
  }
  
  if (d==4)
  {
    if (!all(SignifFeatureData==FALSE))
    {
      x.sig <- x[SignifFeatureData,]
      if (is.vector(x.sig))
        x.sig <- matrix(x.sig, nrow=1)
      
      points3d(x.sig[,1],x.sig[,2], x.sig[,3], size=3,color=featureCol, alpha=trans.alpha)  
    }
  }
}

########## End of addSignifFeatureData ##########
