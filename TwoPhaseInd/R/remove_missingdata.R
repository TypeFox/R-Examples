#!/usr/bin/env Rscript

remove_missingdata <- function(data)
{
  for (i in 1:ncol(data))
  {
    idx <- which(data[,i]=="")
    if (length(idx)>0) data[idx,i] <- NA
  }
  cols <- NULL
  if (sum(is.na(data))>0)
  {
    idx <- NULL
    for (i in 1:ncol(data))
    {
      if (any(is.na(data[,i])))
      {
        idx <- c(idx,which(is.na(data[,i])==TRUE))
        cols <- c(cols,i)
      }
    }
    idx <- unique(idx)
    missingcols <- NULL
    for (i in 1:length(cols))
    {
      missingcols=paste0(missingcols," ",colnames(data)[cols[i]])
    }
    warning(paste0(length(idx), " rows were removed due to missing data in", missingcols))
  }
  
  idxs <- rep(TRUE,nrow(data))
  idxs[idx] <- FALSE
  data <- data[idxs,]
  res <- list()
  res$data <- data
  res$idx <- idxs
  
  return(res)
}
