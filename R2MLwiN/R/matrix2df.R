#' For multiple membership models, translates matrix into a data.frame formatted for MLwiN
#'  
#' Translates a \code{\link[base]{matrix}} into a form usable by MLwiN for multiple membership models,
#' namely a \code{\link[base]{data.frame}} with (a) columns containing membership IDs (if first row matrix is
#' \code{0 1 1 0 1 1}, then first row of generated ID vectors would be, say, \code{2, 3, 5, 6})
#' and (b) columns containing weights (in this example, if \code{standardise = TRUE}, then first
#' row of generated weight vectors would be, say, \code{0.25, 0.25, 0.25, 0.25}, otherwise first
#' row of generated weight vectors would be, say, \code{1, 1, 1, 1}).
#' 
#' @param mat A matrix.
#' @param standardise If \code{TRUE}, ensures the row sums to one; defaults to \code{FALSE}.
#' @param idstub Prefix for columns containing IDs; defaults to \code{id}.
#' @param weightstub Prefix for columns containing weights; defaults to \code{weight}.
#' 
#' @author Zhang, Z., Charlton, C.M.J., Parker, R.M.A., Leckie, G., and Browne,
#' W.J. (2015) Centre for Multilevel Modelling, University of Bristol, U.K.
#' 
#' @seealso \code{\link{df2matrix}}
#' @export
matrix2df <- function(mat, standardise = FALSE, idstub = "id", weightstub = "weight") {
  if (!is.matrix(mat) && !is(mat, "sparseMatrix")) {
    stop("Invalid input data")
  }
  if (is.null(colnames(mat))) {
    colnames(mat) <- 1:ncol(mat)
  }
  if (isTRUE(standardise)) {
    denom <- Matrix::rowSums(mat)
    mat[denom != 0, ] <- mat[denom != 0, ]/denom[denom != 0]
  }
  numvars <- max(Matrix::rowSums(mat != 0))
  idcols <- data.frame(matrix(0, nrow(mat), numvars))
  rownames(idcols) <- 1:nrow(mat)
  colnames(idcols) <- paste0(idstub, 1:numvars)
  weightcols <- data.frame(matrix(0, nrow(mat), numvars))
  colnames(weightcols) <- paste0(weightstub, 1:numvars)
  rownames(weightcols) <- 1:nrow(mat)
  
  for (i in 1:nrow(mat)) {
    row <- mat[i, mat[i, ] != 0, drop = FALSE]
    if (length(row) > 0) {
      idcols[i, 1:length(row)] <- colnames(row)
      weightcols[i, 1:length(row)] <- as.numeric(row)
    }
  }
  
  cbind(idcols, weightcols)
} 
