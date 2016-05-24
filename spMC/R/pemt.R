pemt <-
function(data, coords, mpoints, which.dire, max.dist, tolerance = pi/8,
         rotation = NULL, mle = "avg") {
  # Compute transition probabilities matrices 2D
  # through no ellispoidal interpolation
  #
  #       data vector of data
  #     coords coordinates matrix
  #    mpoints number of points per axes
  # which.dire two choosen 1D directions
  #   max.dist vector of maximum distances
  #  tolerance angle tolerance (in radians)
  #   rotation vector of rotation angles (in radians) to pass to multi_tpfit
  #        mle argument to pass to the function tpfit

  if (!is.factor(data)) data <- as.factor(data)
  if (!is.matrix(coords)) coords <- as.matrix(coords)
  coordsnames <- colnames(coords)
  n <- dim(coords)[1]
  nc <- dim(coords)[2]
  nk <- nlevels(data)

  if (missing(mpoints)) stop("\"mpoints\" is missing")
  if (length(mpoints) != 1) mpoints <- mpoints[1]
  all <- FALSE
  if (missing(which.dire)) all <- TRUE
  if (missing(max.dist)) max.dist <- 1
  if (length(max.dist) == 1) max.dist <- rep(max.dist, nc)
  if (length(max.dist) != nc) 
    stop("\"max.dist\" must be a scalar or vector with the same dimension of coordinates")
  max.dist <- abs(max.dist)

  rawLags <- sapply(max.dist, function(h)
    seq(-h, h, length = mpoints)
  )
  prop <- table(data)
  prop <- prop / sum(prop)

  which.dire <- if (all) combn(1:nc, 2) else matrix(which.dire)
  nimg <- dim(which.dire)
  if (nimg[1] != 2) stop("wrong length of \"which.dire\"")
  if (!is.null(rotation)) {
    if (length(rotation) != nc - 1) stop("wrong length of rotation vector; the right length is ", nc - 1)
  }

  oldcontour <- TRUE
  storage.mode(which.dire) <- "integer"

  imgArgs <- list()
  length(imgArgs) <- nimg[2]
  for (i in 1:nimg[2]) {
    args <- list()
    lagsMat <- as.list(rep(0, nc))
    lagsMat[[which.dire[1, i]]] <- rawLags[, which.dire[1, i]]
    lagsMat[[which.dire[2, i]]] <- rawLags[, which.dire[2, i]]
    lagsMat <- as.matrix(expand.grid(lagsMat))
    nl <- dim(lagsMat)[1]
    
    contour <- oldcontour
    if (contour) {
      x <- multi_tpfit_ml(data, coords, tolerance, rotation)
      if (all(is.finite(unlist(x$coefficients)))) {
        Tprobs <- predict(x, lagsMat)
      }
      else {
        contour <- FALSE
        Tprobs <- NULL
      }
    }

    lagsMat <- t(lagsMat)
    same.dire <- vector("integer", nl)
    same.dire <- .C('wd', lags = as.double(lagsMat), nc = as.integer(nc),
                    nr = as.integer(nl), res = as.integer(same.dire),
                    PACKAGE = "spMC")$res
    wsd <- !duplicated(same.dire)
    coefs <- apply(lagsMat[, same.dire[wsd]], 2,
                   function(d) {
                     res <- tpfit_ml(data, coords, d, tolerance, mle)$coefficients
                     if(prod(is.finite(res)) == 0) res <- matrix(NaN, nk, nk)
                     return(res)
                   })
    storage.mode(wsd) <- "integer"
    nmat <- sum(wsd)
    coefs <- array(unlist(coefs), dim = c(nk, nk, nmat))
    wsd <- as.integer(as.factor(same.dire))
    Eprobs <- array(0, dim = c(nk, nk, nl))
    Eprobs <- .C('predPSEUDO', coefs = as.double(coefs), prop = as.double(prop),
                 lags = as.double(lagsMat), nk = as.integer(nk), nc = as.integer(nc),
                 nr = as.integer(nl), nmat = as.integer(nmat), wsd = as.integer(wsd),
                 whichd = as.integer(which.dire[2, i]), mypred = as.double(Eprobs),
                 NAOK = TRUE, PACKAGE = "spMC")$mypred
    Eprobs <- array(Eprobs, dim = c(nk, nk, nl))

    args$Tprobs <- Tprobs
    args$Eprobs <- Eprobs
    args$X <- rawLags[, which.dire[1, i]]
    args$Y <- rawLags[, which.dire[2, i]]
    imgArgs[[i]] <- args
  }
  imgArgs[[i + 1]] <- list(nk = nk, nomi = levels(data), coordsnames = colnames(coords),
                           mpoints = mpoints, which.dire = which.dire)
  class(imgArgs) <- "pemt"
  return(imgArgs)
}

