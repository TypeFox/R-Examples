getRegions <- function(x) {
  if(sum(substr(dimnames(x)[[1]],4,4)!=".")>0) { #not all regions have 3-character names (need to use slow method)
    output <- unique(as.vector(as.matrix(cbind.data.frame(strsplit(dimnames(x)[[1]],'\\.'))[1,])))
  } else {  #region names all have 3 characters -> fast method
    output <- unique(substr(dimnames(x)[[1]],1,3))
  }
  return(output)
}

"getRegions<-" <- function(x,value) {
  reg <- getRegions(x)
  if(!grepl(".",reg[1],fixed=TRUE)) {
    getCells(x) <- value
    return(x)
  }
  if(length(reg)!=length(value)) stop("Number of regions must agree with current number of regions!")
  tmp <- paste("SAVEREPLACE",dimnames(x)[[1]])
  for(i in 1:nregions(x)) {
    tmp <- sub(paste("SAVEREPLACE ",reg[i],"\\.",sep=""),paste(value[i],"\\.",sep=""),tmp)
  }
  dimnames(x)[[1]] <- tmp
  return(x)
}
