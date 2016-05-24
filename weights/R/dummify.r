dummify <- function(x, show.na=FALSE, keep.na=FALSE){
  if(!is.factor(x)){
    stop("variable needs to be a factor")
  }
  levels(x) <- c(levels(x), "NAFACTOR")
  x[is.na(x)] <- "NAFACTOR"
  levs <- levels(x)
  out <- model.matrix(~x-1)
  colnames(out) <- levs
  attributes(out)$assign <- NULL
  attributes(out)$contrasts <- NULL
  if(show.na==FALSE)
    out <- out[,(colnames(out)=="NAFACTOR")==FALSE]
  if(keep.na==TRUE)
    out[x=="NAFACTOR",] <- matrix(NA, sum(x=="NAFACTOR"), dim(out)[2])
  out
}
