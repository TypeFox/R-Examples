#' Compute Goodman and Kruskal's tau for a dataframe.
#'
#' \code{GKtauDataframe} returns the square matrix of Goodman and Kruskal
#' measures computed between each pair of columns in a dataframe.  Numeric
#' variables in the dataframe are treated as factors.
#'
#' The Goodman and Kruskal tau measure is an asymmetric association measure
#' between two categorical variables, based on the extent to which variation
#' in one variable can be explained by the other.  This function returns an
#' S3 object of class 'GKtauMatrix' that gives the number of levels for
#' each variable on the diagonal of the matrix and the association between
#' variables in the off-diagonal elements.  Note that this matrix is
#' generally NOT symmetric, in contrast to standard correlation matrices.
#'
#' @param df Dataframe from which to compute association measures.
#' @inheritParams GKtau
#' @inheritParams GKtau
#' @return An S3 object of class 'GKtauMatrix' consisting of a square
#' matrix with one row and column for each column of the dataframe df.
#' The structure of this matrix is:
#' \itemize{
#'   \item row and column names are the names of the variables in the dataframe.
#'   \item the diagonal matrix element contains the number of unique levels for
#'   the corresponding variable.
#'   \item off-diagonal matrix elements contain the forward Goodman-Kruskal
#'   tau association from the variable listed in the row names to the
#'   variable listed in the column names.
#' }
#'
#' @author Ron Pearson
#' @export
#'
GKtauDataframe <- function(df, dgts = 3, includeNA = "ifany"){
  #
  n <- ncol(df)
  if (n == 1){
    stop("Dataframe df must have at least two columns")
  }
  GKmatrix <- matrix(nrow = n, ncol = n)
  #
  for (i in 1:n){
    x <- df[, i]
    for (j in 2:n){
      y <- df[, j]
      z <- GKtau(x, y, dgts = dgts, includeNA = includeNA)
      GKmatrix[i, j] <- z$tauxy
      GKmatrix[j, i] <- z$tauyx
    }
    GKmatrix[i, i] <- z$Nx
  }
  colnames(GKmatrix) <- colnames(df)
  rownames(GKmatrix) <- colnames(df)
  class(GKmatrix) <- 'GKtauMatrix'
  return(GKmatrix)
}
