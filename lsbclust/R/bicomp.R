#' Bilinear Decomposition of a Matrix
#' 
#' Decomposes a matrix into an overall mean matrix, row margins matrix, column margins matrix
#' and an interaction matrix, depending on \code{delta}.
#' 
#' @param x A matrix to be decomposed.
#' @param delta A vector of length four with 0/1 entries which controls the type of
#' decomposition made.
#' @param which A vector giving the elements to return, with 0 = original data, 1 = overall means, 
#' 2 = row means, 3 = column means and 4 = interactions.
#' @return An object of class \code{bicomp}, possible also inheriting from class \code{data.frame}, 
#' which is either a named list with the required components, or a single matrix if a single 
#' component is requested. An additional attribute \code{return_type} gives information on the type 
#' of matrices returned.
#' @export
bicomp <- function(x, delta = c(1, 1, 1, 1), which = 0L:4L) {
  
  ## Sanity checks
  delta <- as.numeric(delta)
  if (length(delta) != 4 || !all(delta %in% 0:1))  stop("Argument 'delta' supplied in an unrecognized format.")
  if (!is.matrix(x)) x <- as.matrix(x)
  if (length(which) > 5 || length(which) < 1 || !all(which %in% 0:4)) 
    stop("Argument 'which' contains unrecognized values.")
  
  ## Preliminaries
  J <- nrow(x)
  K <- ncol(x)
  cmatJ.d1 <- if (delta[1]) cmat(J) else diag(J)
  cmatK.d2 <- if (delta[2]) cmat(K) else diag(K)
  cmatK.d3 <- if (delta[3]) cmat(K) else diag(K)
  cmatJ.d4 <- if (delta[4]) cmat(J) else diag(J)
  
  ## Overall mean
  overall <- (delta[1] * delta[3] + delta[2] * delta[4] - delta[1] * delta[2]) * mean(x) * matrix(1, J, K)
  
  ## Row margins
  rows <- delta[2] * cmatJ.d4 %*% x %*% matrix(1, K, K) / K
  
  ## Column margins
  columns <- delta[1] * matrix(1, J, J) %*% x %*% cmatK.d3 / J
  
  ## Interactions
  interactions <- cmatJ.d1 %*% x %*% cmatK.d2
  
  out <- list(original = x, overall = overall, rows = rows, columns = columns, interactions = interactions)
  out <- lapply(out, 'dimnames<-', dimnames(x))
  out <- if (length(which) > 1) out[which + 1] else out[[which + 1]]
  if (length(which) > 1) {
    class(out) <- "bicomp"
  } else class(out) <- c("bicomp", "matrix")
  attr(out, "return_type") <- c("original", "overall", "rows", "columns", "interactions")[which + 1]
  return(out)
}