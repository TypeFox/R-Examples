#' @include utils.R
{}

#for R CMD CHECK
utils::globalVariables(c("nc", "xu"))

#' Generate a reparameterized P-spline base
#'
#' The returned matrix is a low-rank approximation of the original P-spline
#' basis (unless \code{decomposition = "asIs"}), that is projected into the
#' complement of the nullspace of the associated penalty (unless
#' \code{centerBase = FALSE}), i.e. for the default second order difference
#' penalty, the resulting basis cannot reproduce linear or constant functions
#' and parameterizes the "wiggly" part of the influence of \code{x} only. This
#' means that it very rarely makes sense to run a model with \code{sm(x)}
#' without also using \code{\link{lin}(x)} or \code{\link{u}(x)}. The projection
#' improves the separability between the linear and smooth parts of the
#' influence of \code{x} and centers the resulting function estimates s.t
#' \eqn{\sum_i f(x_i) = 0}.
#'
#' @param x covariate
#' @param K number of basis functions in the original basis (defaults to 20)
#' @param spline.degree defaults to 3 for cubic B-plines
#' @param diff.ord order of the difference penalty, defaults to 2 for penalizing
#'   deviations from linearity
#' @param rankZ how many eigenvectors to retain from the eigen decomposition:
#'   either a number > 3 or the proportion of the sum of eigenvalues the
#'   retained eigenvectors must represent at least. Defaults to .999.
#' @param centerBase project the basis of the penalized part into the complement
#'   of the column space of the basis of the unpenalized part? defaults to TRUE
#' @param centerx vector of x-values used for centering (defaults to \code{x})
#' @param decomposition  use a truncated spectral decomposition of the implied
#'   prior covariance of \eqn{f(x)} for a low rank representation with
#'   orthogonal basis functions and i.i.d. coefficients (\code{"ortho"}), or use
#'   the mixed model reparameterization for non-orthogonal basis functions and
#'   i.i.d. coefficients (\code{"MM"}) or use basis functions as they are with
#'   i.i.d. coefficients (\code{"asIs"}). Defaults to \code{"ortho"}.
#' @param tol count eigenvalues smaller than this as zero
#' @author Fabian Scheipl
#' @return a basis for a smooth function in x
#' @references Kneib, T. (2006). Mixed model based inference in structured
#'   additive regression. Dr. Hut.
#'   \url{http://edoc.ub.uni-muenchen.de/archive/00005011/}
#' @importFrom splines  spline.des
#' @importFrom MASS ginv
#' @export
sm <- function(x,
  K = min(length(unique(x)), 20),
  spline.degree = 3, diff.ord = 2, rankZ =.999,
  centerBase = T, centerx = x, decomposition = c("ortho","MM", "asIs"),
  tol = 1e-10) {
  base <- psBasis(x = x, K = K, spline.degree = spline.degree,
    diff.ord = diff.ord)
  if(centerBase) {
    base$X <- centerBase(base$X, x = centerx, degree = diff.ord-1)
  }

  B <- switch(match.arg(decomposition),
    ortho = {
      tmp <- orthoDesign(B = base$X, Cov = ginv(base$P), tol = tol,
        rankZ = rankZ)
      scaleMat(tmp)
    },
    MM = {
      tmp <- mmDesign(B = base$X, K = base$P, tol = tol, rankZ = rankZ)
      scaleMat(tmp)
    },
    asIs = {
      tmp <- base$X
      scaleMat(tmp)
    })
  label <- paste("sm(", deparse(match.call()$x),")", sep ="")
  colnames(B) <- paste(label,".", 1:NCOL(B), sep ="")

  rng <- range(x)
  nd <- !duplicated(x)
  Bu <- B[nd,]
  xu <- x[nd]
  makeNewX <- function(xnew) {
    if(any(xnew<rng[1] | xnew>rng[2])) {
      warning("predictions outside fitted range for ", label,".")
    }
    X <- apply(Bu, 2, function(y) {
      splinefun(x = xu, y = y, 'natural')(xnew)
    })
    return(X)
  }
  makeNewX <- customFunctionEnv(makeNewX, xu = xu, Bu = Bu, rng = rng,
    label = label)
  ret <- structure(B, label = label, predvars = list(makeNewX = makeNewX))
  return(ret)
}
#' Generate orthogonal polynomial base for a numeric covariate without intercept
#'
#' @param x covariate
#' @param order order of the polynomial, defaults to 1
#' @author Fabian Scheipl
#' @return an orthogonal basis for a centered polynomial function in x
#' @export
lin <- function(x, order = 1) {
  B <- poly(x, degree = order)
  label <- paste("lin(", match.call()$x,")", sep ="")
  colnames(B) <- if(order>1) {
    paste(label,".", 1:NCOL(B), sep ="")
  } else {
    label
  }
  scale <- 1/(2 * frob(B))
  B <- scale * B

  rng <- range(x)
  coefs <- attr(B, "coefs")
  nd <- !duplicated(x)
  xu <- x[nd]
  makeNewX <- function(xnew) {
    if(any(xnew < rng[1] | xnew > rng[2])) {
      warning("predictions outside fitted range for ", label,".")
    }
    t(scale * t(poly(xnew, degree = order, coefs = coefs)))
  }
  makeNewX <- customFunctionEnv(makeNewX, coefs = coefs, scale = scale,
    order = order, rng = rng, label = label, xu = xu)

  ret <- structure(B, label = label, predvars = list(makeNewX = makeNewX))
  return(ret)
}

#' Generate design for a factor covariate
#'
#' @param x a covariate that can be treated as a factor
#' @param contr (character) name of the \code{\link[stats]{contrasts}}
#'   associated with the factor. Defaults to \code{"contr.sum"}.
#' @author Fabian Scheipl
#' @return design matrix for a factor
#' @export
fct <- function(x, contr ="contr.sum") {
  x <- factor(x, exclude = NULL)
  contrasts(x) <- do.call(contr, list(n = nlevels(x)))
  B <- model.matrix(~x)[,-1, drop = FALSE]
  scale <- 1/(2 * frob(B))
  B <- scale * B

  label <- paste("fct(", match.call()$x,")", sep ="")
  colnames(B) <- paste(label,".", 1:NCOL(B), sep ="")

  lvls <- levels(x)
  xu <- unique(x)
  makeNewX <- function(xnew) {
    if(!all(xnew %in% lvls)) stop("unknown levels for ", label, ".")
    if(!is.factor(xnew)) xnew <- factor(xnew, levels = lvls)
    contrasts(xnew) <- do.call(contr, list(n = length(lvls)))
    X <- scale * model.matrix(~xnew)[,-1, drop = FALSE]
    return(X)
  }
  makeNewX <- customFunctionEnv(makeNewX, lvls = lvls, contr = contr,
    scale = scale, label = label, xu = xu)
  return(structure(B, label = label, predvars = list(makeNewX = makeNewX)))
}

#' Generate design for a random intercept
#'
#' @param x a covariate that can be treated as a factor
#' @param C an optional known correlation structure of the random effects
#'   associated with x
#' @author Fabian Scheipl
#' @return design matrix for a random intercept for x
#' @export
rnd <- function(x, C = NULL) {
  x <- factor(x, exclude = NULL)
  contrasts(x) <- contr.treatment(n = nlevels(x))
  B <- model.matrix(~ -1 + x)[, , drop = FALSE]
  if(!is.null(C)) {
    cC <- chol(C)
    B <- B %*% t(cC)
  } else cC <- NULL
  scale <- 1/(2 * frob(B))
  B <- scale * B

  label <- paste("rnd(", match.call()$x,")", sep ="")
  colnames(B) <- paste(label,".", 1:NCOL(B), sep ="")

  lvls <- levels(x)
  nd <- !duplicated(x)
  lvlMap <- data.frame(lvl = x[nd], B[nd,])
  makeNewX <- function(xnew) {
    if(!all(xnew %in% lvls)) stop("unknown levels for random effect ", label, ".")
    B <- model.matrix(~ -1 + xnew)[,, drop = FALSE]
    if(!is.null(cC)) {
      B <- B %*% t(cC)
    }
    B <- scale * B
    return(B)
  }
  makeNewX <- customFunctionEnv(makeNewX, lvls = lvls, cC = cC, xu = lvlMap$lvl,
    scale = scale, label = label)
  return(structure(B, label = label, predvars = list(makeNewX = makeNewX)))
}

#' Generate design for an always included covariate
#'
#' Basically a wrapper for \code{model.matrix(~ x, ...)}.
#'
#' @param x covariate
#' @param ... arguments passed to \code{\link{model.matrix}}
#' @author Fabian Scheipl
#' @return a design matrix for x
#' @export
u <- function(x, ...) {
  B <- model.matrix(~ x, ...)[, -1, drop = FALSE]
  label <- paste("u(", match.call()$x,")", sep ="")
  colnames(B) <- if(NCOL(B)>1) {
    paste(label,".", 1:NCOL(B), sep ="")
  } else {
    label
  }
  isFctr <- is.factor(x)
  rng <- if(!isFctr) {
    range(x)
  } else {
    levels(x)
  }

  xu <- unique(x)
  makeNewX <- function(xnew) {
    if(isFctr) {
      if(!all(xnew %in% rng)) stop("unknown levels for ", label, ".")
      xnew <- factor(xnew, levels = rng)
    } else {
      if(any(xnew<rng[1] | xnew>rng[2]))
        warning("predictions outside fitted range for ", label, ".")
    }

    return(model.matrix(~ xnew)[, -1, drop = FALSE])
  }
  makeNewX <- customFunctionEnv(makeNewX, rng = rng, isFctr = isFctr,
    label = label, xu = xu)
  return(structure(B, label = label, predvars = list(makeNewX = makeNewX)))
}
#' Generate design for a 2-D Gaussian Markov Random Field
#'
#' The returned design is (a low-rank approximation to) the matrix square root
#' of the implied covariance of the centered MRF. The function stops if
#' 'islands', i.e. regions without any neighbors are found. Regions without
#' observations have to be removed from the neighborhood matrix and there is
#' currently no \code{predict}-functionality for regions without observations in
#' the original data.
#'
#'
#' @param x a factor: which observation belongs to which region
#' @param N the neighborhood (adjacency) matrix: a symmetric matrix with one
#'   column/row for every level of \code{x}, defining the neighborhood structure
#'   (either 0-1 or with positive weights, e.g. based on shared boundary length
#'   or centroid distances). Has to have rownames and column names that
#'   correspond to the levels of \code{x}, the function checks whether the
#'   rows/columns are in the same order as the levels of \code{x}. Entries on
#'   the diagonal are ignored.
#' @param decomposition  use a (truncated) spectral decomposition of the implied
#'   prior covariance of \eqn{f(x)} for a low rank representation with
#'   orthogonal  basis functions and i.i.d. coefficients (\code{"ortho"}), or
#'   use the mixed model reparameterization for non-orthogonal basis functions
#'   and i.i.d. coefficients (\code{"MM"}). Defaults to \code{"MM"}.
#' @param tol count singular/eigenvalues smaller than this as zero
#' @param rankZ how many eigenvectors to retain from the eigen decomposition:
#'   either a number > 3 or the proportion of the sum of eigenvalues the
#'   retained eigenvectors must represent at least. Defaults to .999.
#' @author Fabian Scheipl
#' @return a transformed design matrix for the Markov Random Field
#' @references Fahrmeir, L., Lang, S. (2001) Bayesian inference for generalized
#'   additive mixed models based on Markov random field priors. \emph{Applied
#'   Statistics}, \bold{50}(2):201--220.
#' @export
mrf <- function(x, N, decomposition = c("ortho", "MM"), tol = 1e-10,
  rankZ =.995) {

  stopifnot(is.factor(x), all(N >= 0), isSymmetric(N), nlevels(x) == ncol(N),
    colnames(N)== rownames(N),
    all(levels(x) == colnames(N)))

  label <- paste("mrf(", match.call()$x,")", sep ="")

  n <- nlevels(x)
  contrasts(x) <- contr.treatment(n = n)


  K <- -N
  diag(K) <- 0
  diag(K) <- -rowSums(K)
  if(any(diag(K)== 0)) {
    warning("region(s) without neighbors found in ", label,
      ", with names: ", paste(colnames(N)[which(diag(K)<1)], collapse =", "),".")
    diag(K)[diag(K)== 0] <- .1
  }

  B <- model.matrix(~ -1 + x)

  B <- switch(match.arg(decomposition),
    ortho = {
      tmp <- orthoDesign(B = B, Cov = ginv(K), tol = tol, rankZ = rankZ)
      scaleMat(tmp)
    },
    MM = {
      tmp <- mmDesign(B = B, K = K, tol = tol, rankZ = rankZ)
      scaleMat(tmp)
    })

  colnames(B) <- paste(label, ".", 1:NCOL(B), sep ="")

  lvls <- levels(x)
  nd <- !duplicated(x)
  lvlMap <- data.frame(lvl = x[nd], B[nd,])
  attributes(lvlMap$lvl)$contrasts <- NULL
  makeNewX <- function(xnew) {
    if(!all(xnew %in% lvls)) stop("unknown regions for ", label,".")
    B <- matrix(0, nrow = length(xnew), ncol = nc)
    for(i in 1:length(lvls)) {
      ind <- which(xnew == lvlMap$lvl[i])
      if(length(ind)) {
        B[ind, ] <- rep(as.numeric(lvlMap[i, -1]), each = length(ind))
      }
    }
    return(B)
  }
  makeNewX <- customFunctionEnv(makeNewX, lvls = lvls, lvlMap = lvlMap,
    xu = lvlMap$lvl, label = label, nc = NCOL(B))
  return(structure(B, label = label, predvars = list(makeNewX = makeNewX)))
}

#' Generate design for penalized surface estimation.
#'
#' This function generates the design for a 2-D penalized spline using (almost)
#' radial basis functions. Use this type of term to account for \emph{spatial}
#' variation. Smooth interactions between covariates are often better fitted via
#' the interactions of \code{\link{lin}} and \code{\link{sm}} terms, because
#' they allow a decomposition into (linear and smooth) marginal trends and
#' (linear-linear, linear-smooth/"varying coefficients", and smooth-smooth)
#' interactions. This decomposition usually makes no sense for spatial data.
#'
#' Note that \code{srf()} expects \code{coords} to be a \code{data.frame} within
#' the larger \code{data.frame} supplied to \code{\link{spikeSlabGAM}} in its
#' \code{data} argument, i.e. \code{coords} is considered a two-dimensional
#' covariate.
#'
#' If \code{baseType} is \code{'thinPlate'}, knot locations for the thin plate
#' spline basis are chosen via a space-filling algorithm (i.e. medoids returned
#' by \code{\link[cluster]{clara}}) as suggested in Ruppert/Wand/Carroll, ch.
#' 13.5. Since the thin plate penalty term penalizes deviations from a linear
#' trend, it is recommended to add marginal linear trends and their interaction
#' to the model if \code{baseType ="thinPlate"} to improve the fit.
#'
#' @note TODO: prediction seems very unstable when extrapolating outside the
#'   convex hull, sometimes even nonsensical values inside conv.hull close to
#'   border - problem in akima::interp?
#'
#' @param coords a \code{data.frame} with two columns containing the coordinates
#' @param K (approximate) number of basis functions in the original basis
#'   (defaults to 50). If \code{baseType ="B"} you can specify a vector giving
#'   the number of marginal basis functions in each direction.
#' @param rankZ how many eigenvectors to retain from the eigen decomposition:
#'   either a number > 3 or the proportion of the sum of eigenvalues the
#'   retained eigenvectors must represent at least. Defaults to .999.
#' @param centerBase project the basis of the penalized part into the complement
#'   of the column space of the basis of the unpenalized part? defaults to TRUE
#' @param baseType Defaults to \code{"B"}, i.e. a tensor product basis based on
#'   marginal cubic B-splines with ridge penalty (i.e. penalizing deviations
#'   from the constant). Set to \code{"thinPlate"} if cubic thin plate splines
#'   are desired, see note below.
#' @param decomposition  use a (truncated) spectral decomposition of the implied
#'   prior covariance of \eqn{f(x, y)} for a low rank representation with
#'   orthogonal	basis functions and i.i.d. coefficients (\code{"ortho"}), or use
#'   the mixed model reparameterization for non-orthogonal basis functions and
#'   i.i.d. coefficients (\code{"MM"}) or use basis functions as they are with
#'   i.i.d. coefficients (\code{"asIs"}). Defaults to \code{"ortho"}.
#' @param tol count eigenvalues smaller than this as zero
#' @author Fabian Scheipl
#' @return a design matrix for the 2-D spline.
#' @references Ruppert, D., Wand, M.P., Carroll, R.J. (2003). Semiparametric
#'   Regression. Cambridge University Press
#' @importFrom cluster clara
#' @importFrom akima interpp
#' @export
srf <- function(coords, K = min(50, sum(nd)/4), rankZ =.999, centerBase = TRUE,
  baseType = c("B", "thinPlate"), decomposition = c("ortho","MM", "asIs"),
  tol = 1e-10) {

  stopifnot(is.data.frame(coords), ncol(coords)== 2)

  nd <- !duplicated(coords)
  baseType <- match.arg(baseType)

  if(baseType == "thinPlate") {
    if(length(K) > 1) {
      warning("Only using first element of K for thin plate basis.")
    }
    base <- tp2DBasis(coords, K[1], nd)
  }
  if(baseType == "B") {
    if(length(K)== 1) K <- c(floor(sqrt(K)), ceiling(sqrt(K)))
    B1 <- psBasis(coords[, 1], K = K[1], spline.degree = 3, diff.ord = 0)$X
    B2 <- psBasis(coords[, 2], K = K[2], spline.degree = 3, diff.ord = 0)$X
    base <- list(X = NULL, P = diag(ncol(B1)* ncol(B2)))
    base$X <- matrix(apply(B1, 2, function(x) {
      x * B2
    }),
      nrow = nrow(B1), ncol = NCOL(B1)* NCOL(B2))
  }
  if(centerBase) {
    base$X <- centerBase(base$X, x = coords[1,], degree = 0)
  }


  B <- switch(match.arg(decomposition),
    ortho = {
      tmp <- orthoDesign(B = base$X, Cov = ginv(base$P), tol = tol,
        rankZ = rankZ)
      scaleMat(tmp)
    },
    MM = {
      tmp <- mmDesign(B = base$X, K = base$P, tol = tol, rankZ = rankZ)
      scaleMat(tmp)
    },
    asIs = base$X)
  label <- paste("srf(", deparse(match.call()$coords),")", sep ="")
  colnames(B) <- paste(label,".", 1:NCOL(B), sep ="")

  rng <- (coords[nd,])[chull(coords[nd,]),]
  Bu <- B[nd,]
  coordsu <- coords[nd,]
  makeNewX <- function(xnew) {
    if(any(!insidePoly(xnew, rng))) {
      warning("Prediction locations outside of convex hull of original observations for ",
        label, ".")
    }

    X <- apply(Bu, 2, function(b) {
      interpp(x = xu[, 1], y = xu[, 2], z = b, xo = xnew[, 1], yo = xnew[, 2],
        extrap = FALSE, linear = FALSE)$z
    })
    return(X)
  }
  makeNewX <- customFunctionEnv(makeNewX, xu = coordsu,  Bu = Bu, rng = rng,
    label = label, insidePoly = insidePoly)
  ret <- structure(B, label = label, rng = rng,
    predvars = list(makeNewX = makeNewX))
  return(ret)
}
