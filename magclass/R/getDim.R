.getDim <- function(elems, x){
  r <- sapply(elems,grepl,fulldim(x)[[2]],fixed=TRUE)
  if(any(colSums(r)==0)) stop("An element was not found in the given data set (",paste(colnames(r)[colSums(r)==0],collapse=", "),")!")
  if(any(colSums(r)>1)) stop("An element was found in more than one dimension in the given data set (",paste(colnames(r)[colSums(r)>1],collapse=", "),"). Please specify the dim to use!")
  if(!any(rowSums(r)==length(elems))) stop("Used elements belong to different dimensions!")
  dim <- which(rowSums(r)==length(elems))
  return(names(fulldim(x)[[2]])[dim])
}