getRegionList <- function(x) {
  return(factor(sub("\\..*$","",dimnames(x)[[1]])))
}

"getRegionList<-" <- function(x,value) {
  reg <- getRegionList(x)
  if(length(reg)!=length(value)) stop("Lengths of RegionLists do not agree!")
  tmp <- sub("^.*\\.","",dimnames(x)[[1]])
  dimnames(x)[[1]] <- paste(as.vector(value),tmp,sep=".")
  return(x)
}
