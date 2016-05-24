cnv <-
function(x, batches, ...){
  if (missing(batches) | NCOL(x)>1)
    return(cnvDefault(x, ...))
  else {
    if (length(unique(batches))==1)
      return(cnvDefault(x, ...))
    else
      return(cnvBatches(x, batches, ...))
  }
}