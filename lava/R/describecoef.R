##' @export
describecoef <- function(x,par,from,to,mean=TRUE) {
  p <- coef(x, mean=mean)
  if (!missing(from)) {
    st1 <- paste0(to,lava.options()$symbol[1],from)
    st2 <- paste0(to,lava.options()$symbol[2],from)
    st3 <- paste0(from,lava.options()$symbol[2],to)
    pos <- na.omit(match(unique(c(st1,st2,st3)),p))
    attributes(pos) <- NULL
    return(pos)
  }
  res <- strsplit(p,lava.options()$symbol[2])
  var.idx <- which(unlist(lapply(res,length))>1) ## Variance parameters
  rest.idx <- setdiff(seq_along(p),var.idx)
  res[rest.idx] <- strsplit(p[rest.idx],lava.options()$symbol[1])
  mean.idx <- which(unlist(lapply(res,length))==1) ## Mean parameters
  reg.idx <- setdiff(rest.idx,mean.idx)
  names(res)[mean.idx] <- paste0("m",seq_along(mean.idx))
  for (i in var.idx)
    attr(res[[i]],"type") <- "cov"
  for (i in mean.idx)
    attr(res[[i]],"type") <- "mean"
  for (i in reg.idx)
    attr(res[[i]],"type") <- "reg"
  if (missing(par))
    return(res)
  return(res[par])
}
