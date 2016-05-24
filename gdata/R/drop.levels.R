
drop.levels  <- function(x, reorder=TRUE, ...)
  UseMethod("drop.levels")

drop.levels.default <- function(x, reorder=TRUE, ...)
  x

drop.levels.factor <- function(x, reorder=TRUE, ...)
{
  x <- x[, drop=TRUE]
  if(reorder) x <- reorder(x, ...)
  x
}

drop.levels.list <- function(x, reorder=TRUE, ...)
{
  lapply(x, drop.levels, reorder=reorder, ...)
}

drop.levels.data.frame <- function(x, reorder=TRUE, ...)
{
  x[] <- drop.levels.list(x, reorder=reorder, ...)
  x
}
