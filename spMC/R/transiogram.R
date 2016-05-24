# transiogram <-
# function(data, coords, direction, max.dist = Inf, mpoints = 20, tolerance = pi/8) {
#   # Empirical transition probabilities estimated by points
#   #
#   #       data vector of data
#   #     coords coordinates matrix
#   #  direction vector (or versor) of choosen direction
#   #   max.dist maximum distance for counting
#   #    mpoints number of lags
#   #  tolerance angle tolerance (in radians)
# 
#   if (!is.numeric(max.dist) || max.dist < 0) stop("\"max.dist\" must be numeric and non negative")
#   if (!is.factor(data)) data <- as.factor(data)
#   if (!is.matrix(coords)) coords <- as.matrix(coords)
#   if (!is.numeric(mpoints) || mpoints <= 0) stop("\"mpoints\" must be positive number")
#   mpoints <- as.integer(ceiling(mpoints))
#   nk <- nlevels(data) # number of categories
#   n <- length(data)   # sample size
#   if (n != dim(coords)[1]) stop("the number of data is not equal to the number of coordinates")
#   nc <- dim(coords)[2]
#   if (length(direction) != nc) stop("wrong length of direction vector")
#   labels <- levels(data)
# 
#   storage.mode(coords) <- "double"
#   storage.mode(direction) <- "double"
# 
#   # definition of bins through I(x <= vDeltaH)
#   direction <- direction / sqrt(sum(direction^2))
#   RNG <- apply(coords, 2, range)
#   dr <- sqrt(sum(diff(RNG)^2))
#   deltah <- min(dr, max.dist) / mpoints
#   vDeltaH <- cumsum(rep(deltah, mpoints))
# 
#   # count transition occurences
#   Tcount <- .C('transCount', n = as.integer(n), data = as.integer(data), 
#                nc = as.integer(nc), coords = as.double(coords),
#                dire = as.double(direction), tolerance = as.double(tolerance), 
#                mpoints = as.integer(mpoints), bins = as.double(vDeltaH), 
#                nk = as.integer(nk), trans = as.double(vector("numeric", nk^2 * mpoints)),
#                DUP = FALSE, PACKAGE = "spMC")$trans
#   Tcount <- array(Tcount, dim = c(nk, nk, mpoints))
#   rwSum <- apply(Tcount, 3, .rowSums, m = nk, n = nk)
#   mtSum <- apply(Tcount, 3, sum)
#   nonComputable <- mtSum == 0
#   if (all(nonComputable)) stop("\"max.dist\" is lower than the minimum distance")
#   # compute transition probabilities
#   Tcount <- .C('transProbs', mpoints = as.integer(mpoints), nk = as.integer(nk), 
#                rwsum = as.double(rwSum), empTR = as.double(Tcount),
#                DUP = FALSE, PACKAGE = "spMC")$empTR
#   res <- list()
#   res$Tmat <- array(Tcount, dim = c(nk, nk, mpoints))
#   res$Tmat <- array(c(diag(, nk), res$Tmat[, , !nonComputable]),
#                     dim = c(nk, nk, mpoints - sum(nonComputable) + 1))
#   colnames(res$Tmat) <- rownames(res$Tmat) <- labels
# 
#   res$lags <- apply(cbind(c(0, vDeltaH[-mpoints]), vDeltaH), 1, mean)[!nonComputable]
#   res$lags <- c(0, res$lags)
#   res$type <- "Empirical"
#   class(res) <- "transiogram"
#   return(res)
# }

transiogram <-
function(data, coords, direction, max.dist = Inf, mpoints = 20, tolerance = pi/8, reverse = FALSE) {
  # Empirical transition probabilities estimated by points
  #
  #       data vector of data
  #     coords coordinates matrix
  #  direction vector (or versor) of choosen direction
  #   max.dist maximum distance for counting
  #    mpoints number of lags
  #  tolerance angle tolerance (in radians)
  #    reverse logical, if TRUE compute probabilities also for the reversible chain

  if (!is.logical(reverse)) stop("\"reverse\" must be a logical value")
  reverse <- reverse[1]
  if (!is.numeric(max.dist) || max.dist < 0) stop("\"max.dist\" must be numeric and non negative")
  if (!is.factor(data)) data <- as.factor(data)
  if (!is.matrix(coords)) coords <- as.matrix(coords)
  if (!is.numeric(mpoints) || mpoints <= 0) stop("\"mpoints\" must be positive number")
  mpoints <- as.integer(ceiling(mpoints))
  nk <- nlevels(data) # number of categories
  n <- length(data)   # sample size
  if (n != dim(coords)[1]) stop("the number of data is not equal to the number of coordinates")
  nc <- dim(coords)[2]
  if (length(direction) != nc) stop("wrong length of direction vector")
  labels <- levels(data)

  storage.mode(coords) <- "double"
  storage.mode(direction) <- "double"

  # definition of bins through I(x <= vDeltaH)
  direction <- direction / sqrt(sum(direction^2))
  RNG <- apply(coords, 2, range)
  dr <- sqrt(sum(diff(RNG)^2))
  deltah <- min(dr, max.dist) / (mpoints * 2^reverse)
  vDeltaH <- cumsum(rep(deltah, mpoints))

  # count transition occurences along the chosen direction
  Tcount <- .C('transCount', n = as.integer(n), data = as.integer(data), 
               nc = as.integer(nc), coords = as.double(coords),
               dire = as.double(direction), tolerance = as.double(tolerance), 
               mpoints = as.integer(mpoints), bins = as.double(vDeltaH), 
               nk = as.integer(nk), trans = as.double(vector("numeric", nk^2 * mpoints)),
               PACKAGE = "spMC")$trans
  Tcount <- array(Tcount, dim = c(nk, nk, mpoints))
  mtSum <- apply(Tcount, 3, sum)
  computable <- mtSum != 0
  computablerev <- c()
  if (reverse) { # count transition occurences along the opposite chosen direction
    revdire <- -1 * direction
    Tcountrev <- vector("numeric", nk^2 * mpoints)
    Tcountrev <- .C('transCount', n = as.integer(n), data = as.integer(data), 
                    nc = as.integer(nc), coords = as.double(coords),
                    dire = as.double(revdire), tolerance = as.double(tolerance), 
                    mpoints = as.integer(mpoints), bins = as.double(vDeltaH), 
                    nk = as.integer(nk), trans = as.double(Tcountrev),
                    PACKAGE = "spMC")$trans
    Tcountrev <- array(Tcountrev, dim = c(nk, nk, mpoints))
    mtSumrev <- apply(Tcountrev, 3, sum)
    computablerev <- mtSumrev != 0
  }
  if (all(!c(computable, computablerev))) stop("\"max.dist\" is lower than the minimum distance")
  # compute transition probabilities
  rwSum <- apply(Tcount, 3, .rowSums, m = nk, n = nk)
  LOSEmat <- .C('transSE', mpoints = as.integer(mpoints), nk = as.integer(nk),
                 rwsum = as.double(rwSum), empTR = as.double(Tcount),
                 se = double(as.integer(mpoints) * nk^2L), PACKAGE = "spMC")$se
  Tcount <- .C('transProbs', mpoints = as.integer(mpoints), nk = as.integer(nk), 
               rwsum = as.double(rwSum), empTR = as.double(Tcount),
               PACKAGE = "spMC")$empTR
  if (reverse) {
    rwSumrev <- apply(Tcountrev, 3, .rowSums, m = nk, n = nk)
    LOSEmatrev <- .C('transSE', mpoints = as.integer(mpoints), nk = as.integer(nk),
                   rwsum = as.double(rwSumrev), empTR = as.double(Tcountrev),
                   se = double(as.integer(mpoints) * nk^2L), PACKAGE = "spMC")$se
    Tcountrev <- .C('transProbs', mpoints = as.integer(mpoints),
                    nk = as.integer(nk), rwsum = as.double(rwSumrev),
                    empTR = as.double(Tcountrev), PACKAGE = "spMC")$empTR
  }
  res <- list()
  if (reverse) {
    wh <- which(c(computablerev, TRUE, computable))
    sq <- seq_len(sum(computablerev))
    wh[sq] <- wh[rev(sq)]
    res$Tmat <- array(c(Tcountrev, diag(, nk), Tcount), dim = c(nk, nk, 1L + 2L * mpoints))[, , wh]
    res$LOSE <- array(c(LOSEmatrev, matrix(0, nk, nk), LOSEmat), dim = c(nk, nk, 1L + 2L * mpoints))[, , wh]
    res$lags <- apply(cbind(c(0, vDeltaH[-mpoints]), vDeltaH), 1, mean)
    res$lags <- c(-rev(res$lags[computablerev]), 0, res$lags[computable])
  }
  else {
    res$Tmat <- array(Tcount, dim = c(nk, nk, mpoints))
    res$Tmat <- array(c(diag(, nk), res$Tmat[, , computable]),
                      dim = c(nk, nk, mpoints - sum(!computable) + 1))
    res$LOSE <- array(LOSEmat, dim = c(nk, nk, mpoints))
    res$LOSE <- array(c(matrix(0, nk, nk), res$LOSE[, , computable]),
                      dim = c(nk, nk, mpoints - sum(!computable) + 1))
    res$lags <- apply(cbind(c(0, vDeltaH[-mpoints]), vDeltaH), 1, mean)[computable]
    res$lags <- c(0, res$lags)
  }
  colnames(res$Tmat) <- rownames(res$Tmat) <- labels

  res$type <- "Empirical"
  class(res) <- "transiogram"
  return(res)
}

