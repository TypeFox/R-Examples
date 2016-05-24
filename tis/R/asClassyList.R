asClassyList <- function(x, ...){
  res <- vector("list", length(x))
  xnames <- names(x)
  names(x) <- NULL
  for (i in seq_along(x)) res[[i]] <- x[i]
  names(res) <- xnames
  res
}
