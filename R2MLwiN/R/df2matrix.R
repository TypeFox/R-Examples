#' Translates a data.frame, formatted for use in multiple membership
#' modelling in MLwiN, to a matrix. 
#' 
#' Translates a \code{\link[base]{data.frame}}, in a form usable by MLwiN for multiple membership models,
#' into a \code{\link[base]{matrix}}. The data.frame needs to contain (a) columns with membership IDs
#' (e.g. first row of which might be \code{2, 3, 5, 6, 0, 0}) and (b) columns containing weights
#' (e.g. first row of which might be \code{0.25, 0.25, 0.25, 0.25, 0, 0}; in this example the first row of
#' resulting matrix would be \code{0, 1, 1, 0, 1, 1}).
#' 
#' @param data A \code{\link[base]{data.frame}} object.
#' @param idcols String vector of the identifier column names.
#' @param weightcols String vector of the weight column names.
#' 
#' @return An adjacency matrix as returned by \code{\link[Matrix]{sparseMatrix}}.
#' 
#' @author Zhang, Z., Charlton, C.M.J., Parker, R.M.A., Leckie, G., and Browne,
#' W.J. (2015) Centre for Multilevel Modelling, University of Bristol, U.K.
#' 
#' @seealso \code{\link{matrix2df}}
#' 
#' @export
df2matrix <- function(data, idcols, weightcols) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("Invalid input data")
  }
  id <- data[, idcols]
  weight <- data[, weightcols]
  id[weight == 0] <- NA
  
  a <- NULL
  for (i in 1:ncol(id)) {
    a <- na.omit(union(a, id[, i]))
  }
  a <- sort(a)
  
  dat <- rep(0, nnzero(weight))
  indi <- dat
  indj <- dat
  
  ind <- 1
  for (i in 1:nrow(id)) {
    for (j in 1:ncol(id)) {
      if (!is.na(id[i, j])) {
        indi[ind] <- i
        indj[ind] <- which(a == id[i, j])
        dat[ind] <- weight[i, j]
        ind <- ind + 1
      }
    }
  }
  
  c <- sparseMatrix(indi, indj, x = dat, dimnames = list(1:nrow(id), a))
  c
} 
