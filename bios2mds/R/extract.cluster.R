extract.cluster <- function (x, align) {

  if (!inherits(x, "kmean"))
    stop("object of class 'kmean' expected")
  if (!inherits(align, "align"))
    stop("object of class 'align' expected")
  nb.clus <- length(x$clusters)
  res <- list ()

  for(i in 1:nb.clus) {
    name <- names(x$clusters[[i]])
    res[[i]] <- align[names(align)%in%name]
    class (res[[i]]) <- c("align")
  }
  names(res)<-sapply(c(1:nb.clus), function(j) {paste("cluster",j, sep = "")})
  return(res)
}