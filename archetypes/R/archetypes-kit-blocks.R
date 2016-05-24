

### Scaling and rescaling functions: #################################

#' Scaling block: standardize to mean 0 and standard deviation 1.
#' @param x Data matrix.
#' @return Standardized data matrix with some attribues.
#' @noRd
std.scalefn <- function(x, ...) {
  m = rowMeans(x)
  x = x - m

  s = apply(x, 1, sd)
  x = x / s

  attr(x, '.Meta') = list(mean=m, sd=s)

  return(x)
}

#' Rescaling block: counterpart of std.scalefn.
#' @param x Standardized data matrix.
#' @param zs Archetypes matrix
#' @return Rescaled archetypes.
#' @noRd
std.rescalefn <- function(x, zs, ...) {

  m = attr(x, '.Meta')$mean
  s = attr(x, '.Meta')$sd

  zs = zs * s
  zs = zs + m

  return(zs)
}



#' Scaling block: no scaling.
#' @param x Data matrix.
#' @return Data matrix.
#' @noRd
no.scalefn <- function(x, ...) {
  return(x)
}

#' Rescaling block: counterpart of no.scalefn.
#' @param x Data matrix.
#' @param zs Archetypes matrix.
#' @return Archetypes zs.
#' @noRd
no.rescalefn <- function(x, zs, ...) {
  if ( is.null(zs) )
    return(matrix(NA, nrow = 0, ncol = 0))

  return(zs)
}



### Dummy and undummy functions: #####################################

#' Dummy block: generator for a dummy function which adds a row
#'   containing a huge value.
#' @param huge The value.
#' @return A function which takes a data matrix and returns the
#'   data matrix with an additonal row containing \code{huge} values.
#' @noRd
make.dummyfn <- function(huge=200) {

  bp.dummyfn <- function(x, ...) {
    y = rbind(x, rep(huge, ncol(x)))

    attr(y, '.Meta') = attr(x, '.Meta')
    attr(y, '.Meta')$dummyrow = nrow(y)

    return(y)
  }

  return(bp.dummyfn)
}


#' Undummy block: remove dummy row.
#' @param x Data matrix.
#' @param zs Archetypes matrix.
#' @return Archetypes zs.
#' @noRd
rm.undummyfn <- function(x, zs, ...) {
  dr = attr(x, '.Meta')$dummyrow

  return(zs[-dr,])
}


#' Dummy block: no dummy row.
#' @param x Data matrix.
#' @return Data matrix x.
#' @noRd
no.dummyfn <- function(x, ...) {
  return(x)
}

#' Undummy block: return archetypes..
#' @param x Data matrix.
#' @param zs Archetypes matrix.
#' @return Archetypes zs.
#' @noRd
no.undummyfn <- function(x, zs, ...) {
  return(zs)
}



### `From X and alpha to archetypes` functions: ######################


#' X to alpha block: QR approach.
#' @param alphas The coefficients.
#' @param x Data matrix.
#' @return The solved linear system.
#' @noRd
qrsolve.zalphasfn <- function(alphas, x, ...) {
  return(t(qr.solve(alphas %*% t(alphas)) %*% alphas %*% t(x)))
}



#' X to alpha block: pseudo-inverse approach.
#' @param alphas The coefficients.
#' @param x Data matrix.
#' @return The solved linear system.
#' @noRd
ginv.zalphasfn <- function(alphas, x, ...) {
  require(MASS)

  return(t(ginv(alphas %*% t(alphas)) %*% alphas %*% t(x)))
}



#' X to alpha block: optim approach.
#' @param alphas The coefficients.
#' @param x Data matrix.
#' @return The solved linear system.
#' @noRd
opt.zalphasfn <- function(alphas, x, ...) {
  z <- rnorm(nrow(x)*nrow(alphas))

  fun <- function(z){
    z <- matrix(z, ncol=nrow(alphas))
    sum( (x - z %*% alphas)^2)
  }

  z <- optim(z, fun, method="BFGS")
  z <- matrix(z$par, ncol=nrow(alphas))

  return(z)
}



### Alpha calculation functions: #####################################


#' Alpha block: plain nnls.
#' @param coefs The coefficients alpha.
#' @param C The archetypes matrix.
#' @param d The data matrix.
#' @return Recalculated alpha.
#' @import nnls
#' @noRd
nnls.alphasfn <- function(coefs, C, d, ...) {
  #require(nnls)

  n = ncol(d)

  for ( j in 1:n )
    coefs[,j] = coef(nnls(C, d[,j]))

  return(coefs)
}

#' Alpha block: nnls with singular value decomposition.
#' @param coefs The coefficients alpha.
#' @param C The archetypes matrix.
#' @param d The data matrix.
#' @return Recalculated alpha.
#' @import nnls
#' @noRd
snnls.alphasfn <- function(coefs, C, d, ...) {
  #require(nnls)

  n = ncol(d)

  nc = ncol(C)
  nr = nrow(C)


  s = svd(C, nv=nc)
  yint = t(s$u) %*% d

  for ( j in 1:n )
    coefs[,j] = coef(nnls(diag(s$d, nrow=nr, ncol=nc) %*% t(s$v), yint[,j]))

  return(coefs)
}



### Beta calculation functions: ######################################


#' Beta block: plain nnls.
#' @param coefs The coefficients beta.
#' @param C The data matrix.
#' @param d The archetypes matrix.
#' @return Recalculated beta.
#' @import nnls
#' @noRd
nnls.betasfn <- nnls.alphasfn



#' Beta block: nnls with singular value decomposition.
#' @param coefs The coefficients beta.
#' @param C The data matrix.
#' @param d The archetypes matrix.
#' @return Recalculated beta.
#' @import nnls
#' @noRd
snnls.betasfn <- snnls.alphasfn



### Norm functions: ##################################################


#' Norm block: standard matrix norm (spectral norm).
#' @param m Matrix.
#' @return The norm.
#' @noRd
norm2.normfn <- function(m, ...) {
  return(max(svd(m)$d))
}


#' Norm block: euclidian norm.
#' @param m Matrix.
#' @return The norm.
#' @noRd
euc.normfn <- function(m, ...) {
  return(sum(apply(m, 2, function(x){sqrt(sum(x^2))})))
}



### Archetypes initialization functions: #############################


#' Init block: generator for random initializtion.
#' @param k The proportion of beta for each archetype.
#' @return A function which returns a list with alpha and beta.
#' @noRd
make.random.initfn <- function(k) {

  bp.initfn <- function(x, p, ...) {

    n <- ncol(x)
    b <- matrix(0, nrow=n, ncol=p)

    for ( i in 1:p )
      b[sample(n, k, replace=FALSE),i] <- 1 / k

    a <- matrix(1, nrow = p, ncol = n) / p

    return(list(betas = b, alphas = a))
  }

  return(bp.initfn)
}

#' Init block: generator for fix initializtion.
#' @param indizes The indizies of data points to use as archetypes.
#' @return A function which returns a list with alpha and beta.
#' @noRd
make.fix.initfn <- function(indizes) {

  fix.initfn <- function(x, p, ...) {
    n <- ncol(x)

    b <- matrix(0, nrow = n, ncol = p)
    b[indizes, ] <- diag(p)

    a <- matrix(1, nrow = p, ncol = n) / p

    return(list(betas = b, alphas = a))
  }

  return(fix.initfn)
}



### Weight functions: ################################################


#' Weight function: move data closer to global center
#' @param data A numeric \eqn{m \times n} data matrix.
#' @param weights Vector of data weights within \eqn{[0, 1]}.
#' @return Weighted data matrix.
#' @noRd
center.weightfn <- function(data, weights, ...) {
  if ( is.null(weights) )
    return(data)

  weights <- as.numeric(1 - weights)

  dr <- attr(data, '.Meta')$dummyrow

  if ( is.null(dr) ) {
    center <- rowMeans(data)
    data <- data + t(weights * t(center - data))
  }
  else {
    center <- rowMeans(data[-dr, ])
    data[-dr, ] <- data[-dr, ] + t(weights * t(center - data[-dr, ]))
  }

  data
}

#' Global weight function: move data closer to global center
#' @param data A numeric \eqn{m \times n} data matrix.
#' @param weights Vector or matrix of data weights within \eqn{[0, 1]}.
#' @return Weighted data matrix.
#' @noRd
center.globweightfn <- function(data, weights, ...) {
  if ( is.null(weights) )
    return(data)

  if ( is.vector(weights) )
    weights <- diag(weights)

  dr <- attr(data, '.Meta')$dummyrow

  if ( is.null(dr) ) {
    data <- data %*% weights
  }
  else {
    data[-dr, ] <- data[-dr, ] %*% weights
  }

  data
}



### Reweights functions: #############################################


#' Reweights function: calculate Bisquare reweights.
#' @param resid A numeric \eqn{m \times n} data matrix.
#' @param reweights Vector of data reweights within \eqn{[0, 1]}.
#' @return Reweights vector.
#' @noRd
bisquare0.reweightsfn <- function(resid, reweights, ...) {
  resid <- apply(resid, 2, function(x) sum(abs(x)))
  resid0 <- resid < sqrt(.Machine$double.eps)

  s <- 6 * median(resid[!resid0])
  v <- resid / s

  ifelse(v < 1, (1 - v^2)^2, 0)
}



#' Reweights function: calculate binary Bisquare reweights.
#' @param resid A numeric \eqn{m \times n} data matrix.
#' @param reweights Vector of data reweights within \eqn{[0, 1]}.
#' @param threshold Threshold for binarization.
#' @return Reweights vector.
#' @noRd
binary.bisquare0.reweightsfn <- function(resid, reweights,
                                         threshold = 0.1, ...) {
  rw <- bisquare0.reweightsfn(resid, reweights, ...)
  ifelse(rw < threshold, 0, 1)
}



### Archetypes family: ###############################################


#' Archetypes family constructor
#'
#' This function returns a problem solving block for each of the
#' different conceptual parts of the algorithm.
#'
#' @param which The kind of archetypes family.
#' @param ... Exchange predefined family blocks with self-defined
#'            functions.
#'
#' @return A list containing a function for each of the different parts.
#'
#' @family archetypes
#'
#' @export
archetypesFamily <- function(which = c('original', 'weighted', 'robust'), ...) {

  which <- match.arg(which)
  blocks <- list(...)

  family <- do.call(sprintf('.%s.archetypesFamily', which), list())
  family$which <- which
  family$which.exchanged <- NULL

  if ( length(blocks) > 0 ) {
    family$which <- sprintf('%s*', family$which)
    family$which.exchanged <- names(blocks)

    for ( n in names(blocks) )
      family[[n]] <- blocks[[n]]
  }


  family
}



#' Original family constructor helper.
#' @return A list of blocks.
#' @noRd
.original.archetypesFamily <- function() {
  list(normfn = norm2.normfn,
       scalefn = std.scalefn,
       rescalefn = std.rescalefn,
       dummyfn = make.dummyfn(200),
       undummyfn = rm.undummyfn,
       initfn = make.random.initfn(1),
       alphasfn = nnls.alphasfn,
       betasfn = nnls.betasfn,
       zalphasfn = qrsolve.zalphasfn,
       globweightfn = function(x, weights) x,
       weightfn = function(x, weights) x,
       reweightsfn = function(x, weights) weights,
       class = NULL)
}

