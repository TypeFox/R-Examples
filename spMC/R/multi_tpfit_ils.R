multi_tpfit_ils <-
function(data, coords, max.dist = Inf, mpoints = 20, tolerance = pi/8, rotation = NULL, q = 10, echo = FALSE, ..., mtpfit) {
  # Estimation for matrixes of transition rates
  #    ( Iterated Least Square Method )
  #
  #       data vector of data
  #     coords coordinates matrix
  #   max.dist maximum distance for counting (expressed for each direction)
  #    mpoints number of lags (expressed for each direction)
  #  tolerance angle for tolerance (in radians, expressed for each direction)
  #   rotation vector of rotation angles (in radians)
  #          q constant greater than one controlling the growth of rho
  #       echo logical value to print the optimization output
  #        ... further arguments to pass to nlminb function
  #     mtpfit multi_tpfit object for a further optimization

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
  if(length(max.dist) < nc) max.dist <- rep(max.dist[1], nc)
  if(length(mpoints) < nc) mpoints <- rep(mpoints[1], nc)
  if(length(tolerance) < nc) tolerance <- rep(tolerance[1], nc)

  res <- list()
  res$coordsnames <- colnames(coords)
  if (is.null(res$coordsnames)) {
    res$coordsnames <- paste("X", 1:nc, sep = "")
  }
  res$coefficients <- list()
  for (i in 1:nc) {
    if (echo) {
      cat("Direction (", sep = "")
      cat(dire.mat[i, ], sep = ", ")
      cat(")\n", sep = "")
    }
    res$coefficients[[i]] <- list()
    if (missing(mtpfit)) {
      res$coefficients[[i]]$coefficients <- tpfit_ils(data, coords, dire.mat[i, ], max.dist[i], 
                                            mpoints[i], tolerance[i], q = q, echo, ...)$coefficients
    }
    else {
      res$coefficients[[i]]$coefficients <- tpfit_ils(data, coords, dire.mat[i, ], max.dist[i], 
          mpoints[i], tolerance[i], q = q, echo, ..., tpfit = mtpfit$coefficients[[i]])$coefficients
    }
    class(res$coefficients[[i]]) <- "tpfit"
    if (echo && i != nc) cat("\n")
  }
  
  res$prop <- table(data)
  res$prop <- res$prop / sum(res$prop)

  res$tolerance <- as.double(tolerance)
  res$rotation <- rotation

  class(res) <- "multi_tpfit"
  return(res)
}
