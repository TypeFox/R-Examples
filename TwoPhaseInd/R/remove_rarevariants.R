#! /usr/bin/env Rscript

remove_rarevariants <- function(data,cutoff=0.02)
{
  #data should be removed if it is a rare variant (freq<cutoff)
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }
  res <- rep(FALSE,ncol(data))
  for (j in 1:ncol(data))
  {
    data1 <- data[,j]
    if (!is.numeric(data1) & ! is.integer(data1)) #not numeric
    {
      if (! is.factor(data1)) data1=as.factor(data1)
      data1=as.numeric(data1)
    }
    counts <- quantile(data1,c(cutoff/2,1-cutoff/2),na.rm=TRUE)
    if (counts[1]==counts[2]) res[j] <- TRUE  #remove
  }
  return(res)
}
