collapseNames <- function(x,collapsedim=NULL) {
  if(is.null(x)) return(NULL)
  f <- fulldim(x)
  if (is.null(collapsedim)) {
    collapsedim <- which(f[[1]][-1:-2]==1)
  } 
  maxdim <- length(f[[1]])-2
  tmp <- getNames(x)
  tmp2 <- names(dimnames(x))[3]
  for(i in collapsedim) {
    searchstring <- paste("^(",paste(rep(".*\\.",i-1),collapse=""),")[^\\.]*(",paste(rep("\\..*",maxdim-i),collapse=""),")$",sep="")
    tmp <- sub(searchstring,"\\1\\2",tmp)
    tmp2 <- sub(searchstring,"\\1\\2",tmp2)
  }
  tmp <- gsub("\\.+","\\.",tmp)
  tmp <- sub("^\\.","",tmp)
  tmp <- sub("\\.$","",tmp)
  tmp2 <- gsub("\\.+","\\.",tmp2)
  tmp2 <- sub("^\\.","",tmp2)
  tmp2 <- sub("\\.$","",tmp2)
  if(length(tmp)==1) if(tmp=="") tmp <- NULL
  if(length(tmp2)==0) tmp2 <- "data"
  getNames(x) <- tmp
  names(dimnames(x))[3] <- tmp2
  x <- clean_magpie(x,what="sets")
  return(x)
}