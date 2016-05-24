#' Hierachical clustering
#' Convenient methods for hierachical clustering
#'
#' @param df data frame
#' @param method method to use, see \code{\link{hclust}}
#' @param metric distance metric to use, see \code{\link{dist}}
#' @param n number of clusters to retrieve, see \code{\link{cut}}
#' @keywords cluster
hierarchical <- function(df, method="complete", metric="euclidean", n=5) {
  if (metric == 'correlation') {
    df <- scale(as.matrix(df))
    metric <- "euclidean"
  }
  as.vector(cutree(hclust(dist(df, metric), method=method), n))
}

#' Clarify matrix
#' Clarify matrix ordering to minimize off diagonals
#'
#' @param a cluster assignments to reassign
#' @param b matrix b
#' @param method clarification method
#' @return vector of reassigned cluster a
#' @keywords manip
#' @seealso \code{\link[e1071]{matchClasses}}
#' @importFrom e1071 matchClasses
clarify <- function(a, b, method="greedy") {
  m <- matchClasses(table(a,b), method=method, verbose=FALSE)
  as.vector(m[a])
}


rescaler <- function(df, type = "sd") {
  f <- switch(type,
    rank = function(x, ...) rank(x, ...),
    var = , sd = function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE),
    robust = function(x) (x - median(x, na.rm = TRUE)) / mad(x, na.rm = TRUE),
    I = function(x) x,
    range = function(x) (x - min(x, na.rm = TRUE)) / diff(range(x, na.rm = TRUE)))

  continuous <- vapply(df, is.numeric, logical(1))
  df[continuous] <- lapply(df[continuous], f)
  df
}

compact <- function(x) Filter(Negate(is.null), x)
