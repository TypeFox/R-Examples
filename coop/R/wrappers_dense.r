is.vec <- function(x)
{
  is.vector(x) && !is.list(x)
}



co_matrix <- function(x, y, type, use)
{
  if (!is.numeric(x))
    stop("argument 'x' must be numeric")
  if (!missing(y))
    stop("argument 'y' can not be used with a matrix 'x'")
  
  use <- check_use(use)
  if (use == "everything")
  {}
  else if (use == "all.obs")
  {
    if (anyNA(x))
      stop("missing observations in covar/pcor/cosine")
  }
  else if (use == "complete.obs")
  {
    if (anyNA(x))
      x <- naomit(x)
  }
  else
    stop("unsupported 'use' method")
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  ret <- .Call(R_co_mat, x, as.integer(type))
  if (!is.null(colnames(x)))
  {
    rownames(ret) <- colnames(x)
    colnames(ret) <- colnames(x)
  }
  
  ret
}



co_vecvec <- function(x, y, type, use)
{
  if (!is.numeric(x))
    stop("argument 'x' must be numeric")
  
  if (missing(y) && type != CO_VAR)
    return(1.0)
  else if (!is.numeric(y))
    stop("argument 'y' must be numeric")
  else if (!is.vec(y))
    stop("argument 'y' must be a non-list vector")
  
  if (length(x) != length(y))
    stop("vectors 'x' and 'y' must have the same length")
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  if (!is.double(y))
    storage.mode(y) <- "double"
  
  # unlike the matrix version, this should come after casting
  # because I don't feel like doing all this garbage for ints
  use <- match.arg(tolower(use), c("everything", "all.obs", "complete.obs"))
  if (use == "everything")
  {}
  else if (use == "all.obs")
  {
    if (anyNA(x) || anyNA(y))
      stop("missing observations in covar/pcor/cosine")
  }
  else if (use == "complete.obs")
  {
    # perhaps a little hacky...
    out <- .Call(R_naomit_vecvec, x, y)
    x <- out[[1]]
    dim(x) <- NULL
    y <- out[[2]]
    dim(y) <- NULL
  }
  
  .Call(R_co_vecvec, x, y, as.integer(type))
}
