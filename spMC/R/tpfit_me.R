tpfit_me <-
function(data, coords, direction, tolerance = pi/8, max.it = 9000, mle = "avg") {
  # Estimation for matrix of transition rates 
  #    ( Maximum Entropy Method )
  #
  #       data vector of data
  #     coords coordinates matrix
  #  direction vector (or versor) of choosen direction
  #     max.it maximum number of iterations for the optimization
  #  tolerance angle tolerance (in radians)
  #        mle argument to pass to the function tpfit

  if (!is.matrix(coords)) coords <- as.matrix(coords)
  n <- dim(coords)[1]
  nc <- dim(coords)[2]
  if (length(direction) != nc) stop("wrong length of direction vector")
  if (!is.factor(data)) data <- as.factor(data)
  nl <- nlevels(data)
  storage.mode(max.it) <- "integer"
  if (n < (nl^2 + nl)) stop("there are not enough data to estimate the parameters")

  proportion <- table(data)
  proportion <- proportion / sum(proportion)

  loc.id <- which_lines(coords, direction, tolerance)
  ml <- mlen(data, coords, loc.id, direction, mle)

  s <- c(proportion) / c(ml)
  fnew <- vector("numeric", nl)
  fnew <- .C('cEmbFrq', s = as.double(s), nk = as.integer(nl), mt = as.integer(max.it),
             eps = as.double(.Machine$double.eps), f = as.double(fnew),
             PACKAGE = "spMC")$f

  Fmat <- outer(fnew, fnew)
  diag(Fmat) <- 0
  Rmat <- diag(1 / (ml * apply(Fmat, 2, sum))) %*%  Fmat
  diag(Rmat) <- -apply(Rmat, 1, sum)

  res <- list()
  res$coefficients <- Rmat
  res$prop <- as.double(proportion)
  names(res$prop) <- levels(data)
  colnames(res$coefficients) <- names(res$prop)
  rownames(res$coefficients) <- names(res$prop)
  res$tolerance <- as.double(tolerance)

  class(res) <- "tpfit"
  return(res)
}

