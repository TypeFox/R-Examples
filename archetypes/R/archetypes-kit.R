#' @include generics.R
#' @include archetypes-kit-blocks.R
#' @include archetypes-class.R
{}



#' Perform archetypal analysis on a data matrix.
#'
#' @param data A numeric \eqn{n \times m} data matrix.
#' @param k The number of archetypes.
#' @param weights Data weights matrix or vector (used as elements of
#'   the diagonal weights matrix).
#' @param maxIterations The maximum number of iterations.
#' @param minImprovement The minimal value of improvement between two
#'   iterations.
#' @param maxKappa The limit of kappa to report an ill-ness warning.
#' @param verbose Print some details during execution.
#' @param saveHistory Save each execution step in an environment for
#'   further analyses.
#' @param family Blocks defining the underlying problem solving mechanisms;
#'   see \code{\link{archetypesFamily}}.
#' @param ... Additional arguments for family blocks.
#'
#' @return An object of class \code{archetypes}, see
#'   \code{\link{as.archetypes}}.
#'
#' @family archetypes
#'
#' @references Cutler and Breiman. Archetypal Analysis. Technometrics,
#'   36(4), 1994. 338-348.
#'
#' @examples
#'   data(toy)
#'   a <- archetypes(toy, 3)
#'
#' @export
archetypes <- function(data, k, weights = NULL, maxIterations = 100,
                       minImprovement = sqrt(.Machine$double.eps),
                       maxKappa = 1000, verbose = FALSE, saveHistory = TRUE,
                       family = archetypesFamily('original'), ...) {

  ### Helpers:
  mycall <- match.call()
  famargs <- list(...)

  memento <- NULL
  snapshot <- function(i) {
    a <- list(archetypes = as.archetypes(t(family$rescalefn(x, family$undummyfn(x, zs))),
              k, alphas = t(alphas), betas = t(betas), rss = rss, kappas = kappas,
              zas = t(family$rescalefn(x, family$undummyfn(x, zas))),
              residuals = resid, reweights = reweights, weights = weights,
              family = list(class = family$class)))

    memento$save(i, a)
  }

  printIter <- function(i) {
    cat(i, ': rss = ', formatC(rss, 8, format = 'f'),
        ', improvement = ', formatC(imp, 8, format = 'f'),
        '\n', sep = '')
  }


  ### Data preparation:
  x1 <- t(data)
  x1 <- family$scalefn(x1, ...)
  x1 <- family$dummyfn(x1, ...)
  x0 <- family$globweightfn(x1, weights, ...)
  x <- x0

  n <- ncol(x)
  m <- nrow(x)


  ### Initialization:
  init <- family$initfn(x, k, ...)

  betas <- init$betas
  alphas <- init$alphas

  zas <- NULL
  zs <- x %*% betas

  resid <- zs %*% alphas - x
  rss <- family$normfn(resid, ...) / n

  reweights <- rep(1, n)

  kappas <- c(alphas=kappa(alphas), betas=kappa(betas),
              zas=-Inf, zs=kappa(zs))
  isIll <- c(kappas) > maxKappa
  errormsg <- NULL

  if ( saveHistory ) {
    memento <- new.memento()
    snapshot(0)
  }


  ### Main loop:
  i <- 1
  imp <- +Inf

  tryCatch(while ( (i <= maxIterations) & (imp >= minImprovement) ) {

    ## Reweight data:
    reweights <- family$reweightsfn(resid, reweights, ...)
    x <- family$weightfn(x0, reweights, ...)


    ## Alpha's:
    alphas <- family$alphasfn(alphas, zs, x, ...)
    zas <- family$zalphasfn(alphas, x, ...)
    rss1 <- family$normfn(zas %*% alphas - x, ...) / n

    kappas[c('alphas', 'zas')] <- c(kappa(alphas), kappa(zas))


    ## Beta's:
    betas <- family$betasfn(betas, x, zas, ...)
    zs <- x %*% betas

    kappas[c('betas', 'zs')] <- c(kappa(betas), kappa(zs))


    ## Residuals, RSS and improvement:
    alphas0 <- family$alphasfn(alphas, zs, x0, ...)

    resid <- zs %*% alphas0 - x0
    rss2 <- family$normfn(resid, ...) / n

    imp <- rss - rss2
    rss <- rss2


    ## Loop Zeugs:
    kappas <- c(alphas = kappa(alphas), betas = kappa(betas),
                zas = kappa(zas), zs = kappa(zs))
    isIll <- isIll & (kappas > maxKappa)

    if ( verbose )
      printIter(i)

    if ( saveHistory )
      snapshot(i)

    i <- i + 1
  },
  error = function(e) errormsg <<- e)


  ### Check illness:
  if ( !is.null(errormsg) ) {
    warning('k=', k, ': ', errormsg)
    return(as.archetypes(NULL, k, NULL, NA, iters = i,
                         call = mycall, history = history,
                         kappas = kappas))
  }

  if ( any(isIll) )
    warning('k=', k, ': ', paste(names(isIll)[isIll], collapse = ', '),
            ' > maxKappa', sep = '')


  ### Rescale and recalculate for original data:
  alphas <- family$alphasfn(alphas, zs, x1)
  betas <- family$betasfn(betas, x1, zs)

  zs <- family$undummyfn(x1, zs)
  zs <- family$rescalefn(x1, zs)

  #resid <- zs %*% alphas - t(data)


  return(as.archetypes(t(zs), k, t(alphas), rss, iters = (i-1),
                       call = mycall, history = memento, kappas = kappas,
                       betas = t(betas), family = family,
                       familyArgs = famargs, residuals = t(resid),
                       weights = weights, reweights = reweights,
                       scaling = attr(x1, ".Meta")))
}



