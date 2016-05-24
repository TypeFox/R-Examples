sample.df <- function (df, size, replace=FALSE, sort=FALSE, prob=NULL) {
  if(!is(df, "data.frame")) stop("first argument must be a data frame")
  N <- nrow(df)
  if (!missing(prob)) {
    if (length(prob) != N) stop("prob= argument must provide weigths for all rows of the data frame")
  }
  idx.sample <- sample.int(N, size, replace=replace, prob=prob)
  if (sort) idx.sample <- sort(idx.sample)
  df[idx.sample, , drop=FALSE]
}
