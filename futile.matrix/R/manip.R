
# Return a portion of a matrix. This is useful for debugging.
peek <- function(x, upper=5, lower=1)
{
  if (is.null(dim(x)))
  {
    my.upper <- min(upper, anylength(x))
    return(x[lower:my.upper])
  }

  upper.row <- min(upper, anylength(x))
  upper.col <- min(upper, ncol(x))
  return(x[lower:upper.row,lower:upper.col])
}


# Expand the matrix m into the larger matrix defined by row.ids x col.ids
# target - either another matrix that defines row and column names or a list
# where target[[0]] are the row names and target[[1]] are the column names
# TODO: Consistency check to ensure all rownames/colnames of m are a subset
# of target
expand <- function(m, target, default=0)
{
  if ('list' %in% class(target)) { dims <- target }
  else if ('matrix' %in% class(target))
  {
    dims <- list(rownames(target), colnames(target))
  }
  else { stop("Argument target must be either a matrix or a list") }

  # Build out rows to new dimensions
  filler.row <- matrix(
    rep(default, ncol(m) * (anylength(dims[[1]]) - nrow(m))), ncol=ncol(m),
    dimnames=list(setdiff(dims[[1]], rownames(m)), colnames(m)) )

  full <- rbind(m, filler.row)
  # Add all other columns
  filler.col <- matrix(
    rep(default, nrow(full) * (anylength(dims[[2]]) - ncol(full))),
    nrow=nrow(full),
    dimnames=list(rownames(full), setdiff(dims[[2]], colnames(full))) )
  full <- cbind(full, filler.col)
  # Sort
  arrange(full)
}

# Select a portion of a matrix based on a regular expression of the row and/or
# column names.
select <- function(m, row.pat=NULL, col.pat=NULL, ...)
{
  out <- m
  if (! is.null(row.pat))
  {
    out <- out[grep(row.pat, rownames(out), ...), , drop=FALSE]
  }

  if (! is.null(col.pat))
  {
    out <- out[, grep(col.pat, colnames(out), ...), drop=FALSE]
  }
  out
}


# Select a portion of a matrix based on a regular expression and assign the
# subset to value. Dimensional integrity is required, otherwise an error will
# result.
"select<-" <- function(m, row.pat=NULL, col.pat=NULL, ..., value)
{
  if (is.null(row.pat) & is.null(col.pat))
  { stop("Either row.pat or col.pat must be set") }

  if (! is.null(row.pat) & is.null(col.pat))
  {
    m[grep(row.pat, rownames(m), ...),] <- value
  }
  else if (is.null(row.pat) & ! is.null(col.pat))
  {
    m[, grep(col.pat, colnames(m), ...)] <- value
  }
  else
  {
    rows <- grep(row.pat, rownames(m), ...)
    cols <- grep(col.pat, colnames(m), ...)
    m[rows,cols] <- value
  }
  invisible(m)
}


# Order rows and columns in a matrix
arrange <- function(m, order.rows=TRUE, order.cols=TRUE, comparator=NULL)
{
  #if (is.null(col.ids)) col.ids <- colnames(m)
  #if (is.null(row.ids)) row.ids <- rownames(m)

  if (order.rows & order.cols)
    m[order(rownames(m)),order(colnames(m)), drop=FALSE]
  else if(order.cols)
    m[,order(colnames(m)), drop=FALSE]
  else if(order.rows)
    m[order(rownames(m)), , drop=FALSE]
}

