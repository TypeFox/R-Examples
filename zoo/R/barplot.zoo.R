barplot.zoo <- function(height, names.arg = NULL, ...)
{
  x <- coredata(height)
  if(!is.null(dim(x))) x <- t(x)
  if(is.null(names.arg)) names.arg <- index2char(index(height))
  barplot(x, names.arg = names.arg, ...)
}
