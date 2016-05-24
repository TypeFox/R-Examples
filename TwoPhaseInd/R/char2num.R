#! /usr/bin/env Rscript


char2num <- function(data)
{
  data <- as.data.frame(data,stringsAsFactors=FALSE)
  res <- data.frame(matrix(NA,nrow=nrow(data),ncol=ncol(data)))
  colnames(res)=colnames(data)
  for (j in 1:ncol(data))
  {
    alllevels <- levels(factor(data[,j]))
    for ( i in 0:(length(alllevels)-1))
    {
      level <- alllevels[i+1]
      idx <- which(data[,j]==level)
      res[idx,j] <- i
    }
  }
  return(res)
}
