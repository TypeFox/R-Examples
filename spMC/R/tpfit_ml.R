tpfit_ml <-
function(data, coords, direction, tolerance = pi/8, mle = "avg") {
  # Estimation for matrix of transition rates
  #          ( Mean Length Method )
  #
  #       data vector of data
  #     coords coordinates matrix
  #  direction vector (or versor) of choosen direction
  #  tolerance angle tolerance (in radians)
  #        mle argument to pass to the function tpfit

  if (!is.factor(data)) data <- as.factor(data)
  if (!is.matrix(coords)) coords <- as.matrix(coords)
  n <- dim(coords)[1]
  nc <- dim(coords)[2]
  if (length(direction) != nc) stop("wrong length of direction vector")
  nl <- nlevels(data)
  if (n < (nl^2 + nl)) stop("there are not enough data to estimate the parameters")

  loc.id <- which_lines(coords, direction, tolerance)
  ml <- mlen(data, coords, loc.id, direction, mle)
  res <- list()
  res$coefficients <- embed_MC(data, coords, loc.id, direction)
  diag(res$coefficients) <- -1
  res$coefficients <- diag(1 / ml) %*% res$coefficients
  res$prop <- table(data)
  res$prop <- as.double(res$prop / sum(res$prop))
  names(res$prop) <- levels(data)
  colnames(res$coefficients) <- names(res$prop)
  rownames(res$coefficients) <- names(res$prop)
  res$tolerance <- as.double(tolerance)
  
  class(res) <- "tpfit"
  return(res)
}
