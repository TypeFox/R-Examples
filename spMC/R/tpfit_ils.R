tpfit_ils <-
function(data, coords, direction, max.dist = Inf, mpoints = 20, tolerance = pi/8, q = 10, echo = FALSE, ..., tpfit) {
  # Estimation for matrix of transition rates 
  #    ( Iterated Least Square Method )
  #
  #       data vector of data
  #     coords coordinates matrix
  #  direction vector (or versor) of choosen direction
  #   max.dist maximum distance for counting
  #    mpoints number of lags
  #  tolerance angle tolerance (in radians)
  #          q constant greater than one controlling the growth of rho
  #       echo logical value to print the optimization output
  #        ... further arguments to pass to nlminb function
  #      tpfit tpfit object for a further optimization

  if (q <= 1) stop("q must be greater than one")
  if (!is.matrix(coords)) coords <- as.matrix(coords)
  n <- dim(coords)[1]
  nc <- dim(coords)[2]
  if (length(direction) != nc) stop("wrong length of direction vector")
  if (!is.factor(data)) data <- as.factor(data)
  nl <- nlevels(data)
  if (n < (nl^2 + nl)) stop("there are not enough data to estimate the parameters")

  EmpTrans <- transiogram(data, coords, direction, max.dist, mpoints, tolerance)
  proportion <- table(data)
  proportion <- proportion / sum(proportion)

  if (missing(tpfit)) {
    Rmat <- matrix(0, nl, nl)
  }
  else {
    Rmat <- tpfit$coefficients
  }
  Bmat <- upper.tri(Rmat) | lower.tri(Rmat)
  Bmat[nl, ] <- FALSE
  lower <- upper <- 0 * Rmat[Bmat]
  upper <- upper + Inf

  optFUN <- function(x, trans, propt, rho = 0) {
    x[is.na(x)] <- 0
    xx <- matrix(0, nl, nl)
    xx[Bmat] <- x
    diag(xx) <- -apply(xx, 1, sum)
    xx[nl,] <- -(propt %*% xx)
    xx[nl,] <- xx[nl,] / propt[nl]
    ll <- length(trans$lags)
    mydim <- c(nl, nl, ll)
    mypred <- array(0, dim = mydim)
    mypred <- .C('predTPFIT', coefficients = as.double(xx),
               prop = as.double(propt), lags = as.double(trans$lags), 
               mydim = as.integer(mydim), mypred = as.double(mypred), 
               PACKAGE = "spMC")$mypred
    sse <- .C('fastrss', n = as.integer(nl * nl * ll), mypred = as.double(mypred),
              Tmat = as.double(trans$Tmat), rss = as.double(0), NAOK = TRUE,
              PACKAGE = "spMC")$rss
    sse <- sse + rho * sum(xx[nl, -nl]^2 * (xx[nl, -nl] < 0), na.rm = TRUE)
    return(sse)
  }

  optmizerFUN <- function(rho, Rmat, Bmat, optFUN, trans, propt, lower, upper, echo, ...) {
    mypar <- as.vector(Rmat[Bmat])
    mypar[is.na(mypar)] <- 0
    res <- nlminb(mypar, optFUN, trans = EmpTrans, propt = c(proportion),
                  rho = rho, lower = as.vector(lower), upper = as.vector(upper),
                  control = list(...))
    if (echo) {
      cat("The algorithm stopped after ")
      cat(res$iterations)
      cat(" steps.\n")
      cat("Penalized sum of square errors: ")
      cat(res$objective)
      cat("\n\n")
    }
    return(res$par)
  }
  
  Rnv <- new.env()
  Rnv <- parent.env(Rnv)
  eps <- sqrt(.Machine$double.eps)
  res <- Rmat[Bmat]

  res <- .Call("bclm", q, eps, res, echo, quote(optmizerFUN(rho, Rmat, Bmat, optFUN, 
        EmpTrans, c(proportion),lower, upper, echo, ...)), Rnv, PACKAGE = "spMC")

  if (echo) cat("Convergence has been reached.\n")

  Rmat <- matrix(0, nl, nl)
  Rmat[Bmat] <- res
  diag(Rmat) <- -apply(Rmat, 1, sum)
  Rmat[nl,] <- -(c(proportion) %*% Rmat)
  Rmat[nl,] <- Rmat[nl,] / proportion[nl]

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
