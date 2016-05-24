multi_tpfit_me <-
function(data, coords, tolerance = pi/8, max.it = 9000, rotation = NULL, mle = "avg") {
  # Estimation for matrixes of transition rates
  #    ( Maximum Entropy Method )
  #
  #       data vector of data
  #     coords coordinates matrix
  #  tolerance angle for tolerance (in radians)
  #     max.it maximum number of iterations for the optimization
  #   rotation vector of rotation angles (in radians)
  #        mle argument to pass to the function tpfit

  if (!is.factor(data)) data <- as.factor(data)
  if (!is.matrix(coords)) coords <- as.matrix(coords)
  n <- dim(coords)[1]
  nc <- dim(coords)[2]
  nk <- nlevels(data)
  if (n < (nk + nk^2 * nc)) stop("there are not enough data to estimate the parameters")
  dire.mat <- diag(, nc)
  if (!is.null(rotation)) {
    if (length(rotation) != nc - 1) stop("wrong length of rotation vector; the right length is ", nc - 1)
    storage.mode(rotation) <- "double"
    dire.mat <- .C('rotaxes', nc = as.integer(nc), ang = as.double(rotation),
                   res = as.double(dire.mat), PACKAGE = "spMC")$res
    dire.mat <- t(matrix(dire.mat, nc, nc))
  }

  res <- list()
  res$coordsnames <- colnames(coords)
  if (is.null(res$coordsnames)) {
    res$coordsnames <- paste("X", 1:nc, sep = "")
  }
  res$coefficients <- apply(dire.mat, 1, function(d) {
    md.coef <- list()
    md.coef$coefficients <- tpfit_me(data, coords, d, tolerance, max.it, mle)$coefficients
    class(md.coef) <- "tpfit"
    return(md.coef)
  })
  
  res$prop <- table(data)
  nameK <- names(res$prop)
  res$prop <- as.double(res$prop / sum(res$prop))
  names(res$prop) <- nameK

  res$tolerance <- as.double(tolerance)
  res$rotation <- rotation

  class(res) <- "multi_tpfit"
  return(res)
}
