long2mitml.list <- function(x, split, exclude=NULL){
# convert data set in "long" format to mitml.list

  i1 <-  which(colnames(x)==split)
  f <- x[,i1]
  if(!is.null(exclude)){
    i2 <- if(length(exclude)==1) f!=exclude else !f%in%exclude
    x <- x[i2,,drop=F]
    f <- f[i2]
  }
  out <- split(x[,-i1,drop=F], f=f)
  names(out) <- NULL

  class(out) <- c("mitml.list","list")
  out

}

