#' A package for hierarchical Bayesian small area estimation.
#'
#' Package hbsae provides functions to compute small area estimates based on the
#' basic unit-level and area-level models. The models are fit and small area estimates are
#' computed in a hierarchical Bayesian way, using numerical integration.
#'
#' The main function that does most of the computational work is \code{\link{fSAE.Unit}}. Function \code{\link{fSAE}} is provided as
#' a more convenient interface to \code{\link{fSurvReg}}, \code{\link{fSAE.Area}} and \code{\link{fSAE.Unit}}.
#'
#' @name hbsae-package
#' @aliases hbsae
#' @docType package
#' @author Harm Jan Boonstra \email{hjboonstra@@gmail.com}
#' @import Matrix
NULL

#' Fit a linear model with random area effects and compute small area estimates.
#'
#' This function prepares the (unit-level) input data and calls one of the lower level functions \code{\link{fSurvReg}}, \code{\link{fSAE.Area}} or \code{\link{fSAE.Unit}}
#' to compute survey regression, area-level model or unit-level model small area estimates. Area-level model estimates 
#' are computed by first computing survey regression estimates and using these as input for \code{\link{fSAE.Area}}.
#'
#' @export
#' @param formula model formula, indicating response variable and covariates.
#' @param data unit-level data frame containing all variables used in \code{formula}, \code{area} and \code{formula.area} arguments.
#'        These variables should not contain missing values.
#' @param area name of area indicator variable in \code{data}; if \code{NULL}, no random effects are used in the model.
#' @param popdata data frame or matrix containing area population totals for all covariates. The rows should correspond to areas
#'        for which estimates are required.
#'        Column names should include those produced by \code{model.matrix(formula, data, contrasts.arg)}, up to permutations of the names in interactions.
#'        A column named '(Intercept)' is required and should contain the area population sizes.
#'        If \code{popdata} is \code{NULL}, only the model fit will be returned.
#' @param type type of small area estimates: "direct" for survey regression, "area" for area-level model, "unit" for unit-level model estimates.
#'        If \code{type} is "data" then only the data including the model matrix and population means are returned.
#' @param model.direct if type="area", this argument can be used to specify by means of a formula the covariates to use for the computation of the initial survey regression estimates.
#'        If unspecified, the covariates specified by \code{formula} are used both at the unit-level (for the initial estimates) and at the area-level (for the area-level model estimates).
#' @param formula.area if type="unit", this is an optional formula specifying covariates that should be used at the area level.
#'        These covariates should be available in \code{popdata}.
#' @param contrasts.arg list for specification of contrasts for factor variables. Passed to \code{model.matrix}.
#' @param remove.redundant if \code{TRUE} redundant columns in the design matrix are removed. A warning is issued if
#'        the same redundancy does not show also in the corresponding population totals. In the case of the
#'        area level model there may still be redundancy at the area level.
#' @param redundancy.tol tolerance for detecting linear dependencies among the columns of the design matrix. Also used as tolerance in the check whether the design matrix reduncancy is shared by the population totals.
#' @param sparse if \code{TRUE} \code{sparse.model.matrix} (package \code{Matrix}) is used to compute the covariate design matrix. This can be efficient
#'        for large datasets and a model containing categorical variables with many categories.
#' @param ... additional arguments passed to \code{\link{fSAE.Unit}} or \code{\link{fSurvReg}}.
#' @return An object of class \code{sae} containing the small area estimates, their MSEs, and the model fit. If \code{type} is "data" a list containing
#'         the model matrix, response vector, area indicator, area population sizes and matrix of population means is returned.
#' @seealso \code{\link{sae-class}}
#' @example examples/fSAE.R
fSAE <- function(formula, data, area=NULL, popdata=NULL, type="unit", model.direct=NULL, formula.area=NULL,
                 contrasts.arg=NULL, remove.redundant=TRUE, redundancy.tol=1e-7, sparse=FALSE, ...) {

  if (is.null(area))
    areavar <- factor(rep.int("global_", nrow(data)))
  else {
    if (!( is.character(area) && (length(area) == 1) )) stop("area must be a single variable name")
    areavar <- na.fail(as.factor(data[[area]]))
    if (!length(areavar)) stop("area variable not found in data")
  }

  type <- casefold(type)
  if (!(type %in% c("direct", "area", "unit", "data"))) stop("type must be one of 'direct', 'area', 'unit' or 'data'")
  if (!is.null(model.direct) && (type != "area")) warning("argument model.direct ignored; only used for area-level model")
  if (!is.null(popdata)) {
    if (is.vector(popdata)) {
      if (nlevels(areavar) == 1)  # a single area
        popdata <- matrix(popdata, nrow=1, dimnames=list("global_", names(popdata)))
	  else  # assume that the vector contains area population totals
        popdata <- matrix(popdata, ncol=1, dimnames=list(names(popdata), "(Intercept)"))
    }
    if (!is.element("(Intercept)", colnames(popdata)))
      stop("popdata must contain a column named '(Intercept)' with area population sizes")
    Narea <- popdata[, "(Intercept)"]  # area population sizes
    # if area is a column in popdata than use it as the area label, otherwise use rownames of popdata
    if (!is.null(area) && (area %in% colnames(popdata))) {
      dimnames(popdata)[[1]] <- names(Narea) <- as.character(as.factor(popdata[, area]))
    } else
      names(Narea) <- rownames(popdata)
  } else {
    Narea <- NULL
  }

  mf <- model.frame(formula, data, na.action=na.fail)
  y <- as.vector(model.response(mf, "numeric"))
  if (is.null(y)) stop("no response variable")
  if (sparse) {
    X <- sparse.model.matrix(formula, mf, contrasts.arg)
  } else {
    # faster, but may consume more memory during computation
    X <- Matrix(model.matrix(formula, mf, contrasts.arg))
  }
  
  if (!is.null(popdata)) {
    # order the interactions alphabetically to match X and popdata
    ma <- match(orderInteractions(colnames(X)), orderInteractions(colnames(popdata)))
    if (any(is.na(ma))) {
      cat("Some covariates could not be found in argument popdata:\n")
      print(colnames(X)[which(is.na(ma))])
      stop("missing columns in popdata")
    }
    PopMeans <- as.matrix(popdata[, ma, drop=FALSE]) / Narea
  } else {
    if (!(type %in% c("unit", "data"))) stop("popdata required for type 'direct' or 'area'")
    PopMeans <- NULL
  }

  if (!is.null(formula.area)) {  # add area-level covariates to X
    if (is.null(popdata) || !all(levels(areavar) %in% rownames(popdata))) stop("popdata for all sampled areas required to derive area level covariates")
    covnames.area <- colnames(model.matrix(formula.area, data[1, ]))
    ma <- match(orderInteractions(covnames.area), orderInteractions(colnames(popdata)))
    if (any(is.na(ma))) {
      cat("Some area-level covariates could not be found in argument popdata:\n")
      print(covnames.area[which(is.na(ma))])
      stop("missing columns in popdata")
    }
    PopMeans.A <- as.matrix(popdata[, ma, drop=FALSE]) / Narea
    dimnames(PopMeans.A)[[2]] <- paste(colnames(PopMeans.A), "A", sep=".")
    PopMeans <- cbind(PopMeans, PopMeans.A)
    X <- cbind2(X, PopMeans.A[as.character(areavar), ])
  }

  columns.removed <- NULL
  if (remove.redundant) {
    # check for redundant covariates
    XX <- as.matrix(crossprod(X))  # for the following check, use ordinary matrix type
    test <- qr(XX, tol=redundancy.tol)  # use qr on XX instead of X for (memory-)performance reasons
    if (test$rank < ncol(X)) {
      rcol <- test$pivot[(test$rank + 1):ncol(X)]
      if (!is.null(popdata)) {
        # check that PopMeans displays same redundancy as X; disable next three lines to skip this test
        temp <- eigen(XX)$vectors[, (test$rank + 1):ncol(X)]
        if (sum(abs(PopMeans %*% temp) > redundancy.tol))
          warning("redundancy in covariate matrix not shared by matrix of population means")
        PopMeans <- PopMeans[, -rcol, drop=FALSE]  # remove redundant columns from PopMeans
      }
      columns.removed <- colnames(X)[rcol]
      X <- X[, -rcol, drop=FALSE]  # remove redundant columns from X
    }
  }

  out <- switch(type,
    direct = fSurvReg(y, X, areavar, Narea, PopMeans, ...),
    area = {
      if (!is.null(model.direct))
        sreg <- fSAE(model.direct, data, area, popdata, type="direct", ...)
      else
        sreg <- fSurvReg(y, X, areavar, Narea, PopMeans, ...)
      fSAE.Area(est.init=sreg$est, var.init=sreg$mse, X=PopMeans, Narea=Narea, ...)
    },
    unit = fSAE.Unit(y, X, areavar, Narea, PopMeans, ...),
    list(y=y, X=X, area=areavar, Narea=Narea, PopMeans=PopMeans)
  )

  out$call <- match.call()
  out$columns.removed <- columns.removed
  out
}


#' Compute small area estimates based on the survey regression estimator.
#'
#' This function computes survey regression estimates as a special case of unit-level model small area estimates with a (relatively) very large value for the between-area variance
#' but without including area effects in the model fit. The model assumes a single overall variance parameter, so that the resulting estimated variances are not area-specific but smoothed.
#' Varying inclusion probabilities may be taken into account by including them in the model, e.g. as an additional covariate,
#' and/or as model variance structure by specifying arguments v and vpop, see \code{\link{fSAE.Unit}}. The resulting estimates may be used as input estimates for area-level model small area estimation.
#'
#' @export
#' @param y response vector of length n.
#' @param X n x p model matrix.
#' @param area n-vector of area codes, typically a factor variable with m levels, where m is the number of sampled areas.
#' @param Narea M-vector of area population sizes.
#' @param Xpop M x p matrix of population means.
#' @param removeEmpty whether out-of-sample areas should be removed from the results. If \code{FALSE} these areas are retained in the vectors of estimates,
#'   but they will have (relatively) very large standard errors.
#' @param ... optional arguments v and vpop passed to \code{\link{fSAE.Unit}}.
#' @return An object of class \code{sae} containing the survey regression small area estimates and their estimated variances.
#' @seealso \code{\link{sae-class}}
#' @references 
#'   G.E. Battese, R.M. Harter and W.A. Fuller (1988).
#'     An Error-Components Model for Prediction of County Crop Areas Using Survey and Satellite Data.
#'     Journal of the American Statistical Association, 83(401), 28-36.
#'
#'   J.N.K. Rao (2003). Small Area Estimation. Wiley.
#' @example examples/fSurvReg.R
fSurvReg <- function(y, X, area, Narea, Xpop, removeEmpty=TRUE, ...) {
  funArgs <- list(...)
  funArgs <- funArgs[na.omit(match(c("v", "vpop", "keep.data"), names(funArgs)))]  # drop unused arguments
  funArgs <- c(list(y=y, X=X, area=area, Narea=Narea, Xpop=Xpop, method="survreg", CV=FALSE, full.cov=FALSE), funArgs)
  out <- do.call(fSAE.Unit, funArgs)
  if (removeEmpty) {  # remove any out-of-sample areas
    compNames <- c("est", "synth", "g1", "g2", "mse", "Narea", "predAreaNames", "wpop")
    out[compNames] <- lapply(out[compNames], function(comp) comp[out$inSampleAreas])
    out$inSampleAreas <- order(out$inSampleAreas)
    out$M <- out$m
  }
  out$type <- "direct"
  out
}


#' Compute small area estimates based on the basic area-level model.
#'
#' This function returns small area estimates based on the basic area-level model, also known as the Fay-Herriot model.
#' It calls \code{\link{fSAE.Unit}} to carry out the computations.
#'
#' @export
#' @param est.init m-vector of initial estimates, where m is the number of sampled areas.
#' @param var.init m-vector of corresponding variance estimates.
#' @param X M x p matrix of area population means, where M is the number of areas for which estimates are required.
#' @param x defaults to X, but in case auxiliary information is available at the unit level, it may be set to the corresponding matrix of covariate sample means.
#' @param ... additional arguments passed to \code{fSAE.Unit}.
#' @return An object of class \code{sae} containing the small area estimates and MSEs, the model fit, and model selection measures.
#' @seealso \code{\link{sae-class}}
#' @references
#'   R.E. Fay and R.A. Herriot (1979). Estimates of Income for Small Places: An Application of James-Stein Procedures to Census Data.
#'   Journal of the American Statistical Association 74(366), 269-277.
#'
#'   J.N.K. Rao (2003). Small Area Estimation. Wiley.
#' @example examples/fSAE.Area.R
fSAE.Area <- function(est.init, var.init, X, x=X, ...) {
  if (nrow(X) > length(est.init))
    x <- x[names(est.init), , drop=FALSE]
  funArgs <- list(...)
  funArgs$v <- funArgs$vpop <- funArgs$nu0 <- funArgs$s20 <- NULL  # these are defined below
  # NB fpc must be FALSE for area-level model
  funArgs <- c(list(y=est.init, X=x, area=if (!is.null(names(est.init))) factor(names(est.init), names(est.init)) else 1:length(est.init),
	                Xpop=X, fpc=FALSE, v=var.init, nu0=1e4*length(est.init), s20=1), funArgs)
  out <- do.call(fSAE.Unit, funArgs)
  out$type <- "area"
  out
}


#' Compute small area estimates based on the basic unit-level model.
#'
#' This is the function that carries out most of the computational work. It computes small area estimates based on the basic unit-level model, also known as the
#' Battese-Harter-Fuller model, although it is also called by \code{\link{fSurvReg}} and \code{\link{fSAE.Area}} to compute survey regression
#' or area-level model small area estimates. By default, Hierarchical Bayes estimates are computed, using fast one-dimensional
#' numerical integration to average over the posterior density for the ratio of between and within area variance. This way, the small area estimates
#' and MSEs account for the uncertainty about this parameter. Besides hierarchical Bayes, REML and hybrid methods are supported.
#' These methods use the REML estimate or posterior mean of the variance ratio, respectively, as a plug-in estimate. Both methods do not account for uncertainty about this
#' parameter. Synthetic estimates are computed by setting the variance ratio to zero.
#'
#' The default Hierarchical Bayes method uses numerical integration (as provided by function \code{\link{integrate}}) to compute
#' small area estimates and MSEs. The model parameters returned, such as fixed and random effects, are currently not averaged over the
#' posterior distribution for the variance ratio. They are evaluated at the posterior mean of the variance ratio.
#'
#' @export
#' @param y response vector of length n.
#' @param X n x p model matrix.
#' @param area n-vector of area codes, typically a factor variable with m levels, where m is the number of sampled areas.
#' @param Narea M-vector of area population sizes, where M is the number of areas for which estimates are required.
#'        There should be a ono-to-one correspondence with the rows of \code{Xpop}.
#'        This argument is required unless \code{Xpop=NULL} or \code{fpc=FALSE}.
#' @param Xpop M x p matrix of population means. If \code{Xpop} is not provided, only the model fit is returned.
#' @param fpc whether a finite population correction should be used. Default is \code{TRUE}.
#' @param v unit-level variance structure, n-vector. Defaults to a vector of 1s. In some cases it might be useful
#'        to take v proportional to the sampling probabilities.
#' @param vpop population area means of v, M-vector. Defaults to a vector of 1s. Not used when \code{fpc} is \code{FALSE}.
#' @param w area-level variance structure, m-vector. Defaults to a vector of 1s.
#' @param wpop area-level variance structure, M-vector. Defaults to a vector of 1s.
#'        Only components of \code{wpop} corresponding to out-of-sample areas are actually used.
#' @param method one of "HB", "hybrid", "REML", "synthetic", "survreg", "BLUP" where
#'        "HB" (default) does the full hierarchical Bayes computation, i.e. numerical integration over the posterior density for the between area variance parameter,
#'        "hybrid" computes the Best Linear Unbiased Predictor (BLUP) with the posterior mean for the variance parameter plugged in,
#'        "REML" computes the BLUP with the restricted maximum likelihood estimate of the variance parameter plugged in,
#'        "synthetic" computes synthetic estimates where the between area variance is set to 0, and
#'        "survreg" computes survey regression estimates where the between area variance approaches infinity.
#'        "BLUP" computes BLUP estimates with the value provided for \code{lambda0} as a fixed plug-in value for the ratio of between and within area variance.
#'        Only method "HB" takes uncertainty about the between-area variance into account.
#' @param beta0 mean vector of normal prior for coefficient vector.
#' @param Omega0 inverse covariance matrix of normal prior for coefficient vector. Default prior
#'        corresponds to the (improper) uniform distribution.
#' @param nu0 degrees of freedom parameter for inverse gamma prior for residual (within-area) variance. Default is 0.
#' @param s20 scale parameter for inverse gamma prior for residual (within-area) variance. Default is 0.
#         The default prior is p(se2) ~ 1/se2.
#' @param prior prior density for the ratio lambda = between-area-variance / within-area variance. This should be a (vectorized) function that takes
#'        a vector lambda and returns a vector of prior density values at lambda. The density does not have to be normalized. The default is
#'        the (improper) uniform prior. The within-area variance and lambda are assumed independent a priori.
#' @param CV whether leave-one-out cross-validation and other model selection measures should be computed. Default \code{TRUE}.
#' @param CVweights n-vector of weights to use for CV computation.
#' @param silent if \code{FALSE}, plot the posterior density for the variance ratio.
#' @param keep.data if \code{TRUE} return the input data (y,X,area,Xpop).
#'        This is required input for the cross-validation function \code{CVArea}.
#' @param full.cov if \code{TRUE} compute the full covariance matrix for the small area estimates.
#'        The computed correlations do not account for uncertainty about the variance ratio.
#' @param lambda0 optional starting value for the ratio of between and within-area variance used in the numerical routines.
#'        If \code{method="BLUP"} then this value will instead be used as a fixed plug-in value.
#' @param rel.int.tol tolerance for the estimated relative integration error (default is 1 percent).
#'        A warning is issued if the estimated relative error exceeds this value.
#' @param ... additional control parameters passed to function \code{integrate}.
#' @return An object of class \code{sae} containing the small area estimates and MSEs, the model fit, and model selection measures.
#' @seealso \code{\link{sae-class}}
#' @references
#'   G.E. Battese, R.M. Harter and W.A. Fuller (1988).
#'     An Error-Components Model for Prediction of County Crop Areas Using Survey and Satellite Data.
#'     Journal of the American Statistical Association, 83(401), 28-36.
#'
#'   G.S. Datta and M. Ghosh (1991). Bayesian Prediction in Linear Models: Applications to Small Area Estimation.
#'     The Annals of Statistics 19(4), 1748-1770.
#'
#'   J.N.K. Rao (2003). Small Area Estimation. Wiley.
#' @example examples/fSAE.Unit.R
fSAE.Unit <- function(y, X, area, Narea=NULL, Xpop=NULL, fpc=TRUE, v=NULL, vpop=NULL, w=NULL, wpop=NULL, method="HB",
                      beta0=rep(0, ncol(X)), Omega0=Diagonal(n=ncol(X), x=0), nu0=0, s20=0, prior=function(x) rep.int(1L, length(x)),
                      CV=TRUE, CVweights=NULL, silent=FALSE, keep.data=FALSE, full.cov = nrow(Xpop) < 1000,
                      lambda0=NULL, rel.int.tol=0.01, ...) {

  if (!(method %in% c("HB", "hybrid", "REML", "synthetic", "survreg", "BLUP"))) stop(sprintf("unrecognized method '%s'", method))
  if (method == "BLUP" & is.null(lambda0)) stop("for method \"BLUP\" a value for lambda0 must be provided")
  
  area <- as.factor(area)  # make sure area is a factor variable
  area <- area[, drop=TRUE]  # drop levels corresponding to unobserved areas
  sampledAreaNames <- levels(area)  # names of sampled areas
  area.aggr <- aggrMatrix(area)  # sparse aggregation matrix for fast aggregation by area
  vn <- as.vector(table(area))  # area sample sizes
  #vn <- as.integer(colSums(area.aggr))  # colSums for ddiMatrix does not work in R <= 2.10
  m <- length(vn)  # number of areas in sample
  n <- length(y)
  
  if (!is.null(Xpop)) {
    M <- nrow(Xpop)  # number of areas for which estimates are computed
    if (is.null(rownames(Xpop)))
      if (M == m) {  # assume same order as levels of area
        warning("no row names in Xpop; assuming same order as levels of area variable")
        rownames(Xpop) <- sampledAreaNames
      } else
        stop("cannot identify row names of Xpop as area names; please add row names to Xpop")
    else
      if ((M >= m) && !all(sampledAreaNames %in% rownames(Xpop)))
        warning("not all sampled area names can be matched to row names of Xpop")
  } else
    M <- 0L

  if (is.null(v)) {
    v <- 1L
    vpop <- rep.int(1L, M)  # areas for which estimates are computed
  }
  # distinguish the equal variances case for optimization purposes
  if (all(v == v[1])) {
    v <- v[1]
    equalv <- TRUE
  } else
    equalv <- FALSE

  if (is.null(w)) {
    w <- rep.int(1L, m)  # observed areas
    wpop <- rep.int(1L, M)  # areas for which estimates are computed
  }

  if (!is(X, "Matrix")) {
    if (is.data.frame(X))
	  X <- as.matrix(X)
    X <- Matrix(X)
  }
  if (!is.null(rownames(X)))
    dimnames(X) <- list(NULL, dimnames(X)[[2]])  # strip row labels, if any
  if (M > 0) Xpop <- as.matrix(Xpop)
  # in the code below, assume X has class Matrix and Xpop class matrix

  # some data integrity checks
  if (nrow(X) != n) stop("incompatible dimensions of y and X")
  if (length(area) != n) stop("arguments area and y have different lengths")
  if (((length(v) != 1) && (length(v) != n)) || (any(v <= 0))) stop("v must be a positive vector of length 1 or length(y)")
  if ((m == 1) && (method %in% c("REML", "hybrid", "HB", "BLUP"))) stop("model without area effects; use method synthetic or survreg")
  if ((length(w) != m) || (any(w <= 0))) stop("w must be a positive vector with length equal to the number of sampled areas")
  if (M > 0) {
    inSampleAreas <- as.vector(na.omit(match(sampledAreaNames, rownames(Xpop))))  # indices of sampled areas in population (prediction) objects
    predSampledAreas <- sampledAreaNames %in% rownames(Xpop)  # sampled areas for which predictions are made (typically all)
    if (ncol(Xpop) != ncol(X)) stop("number of columns of Xpop and X differ")
    if ((length(wpop) != M) || (any(wpop <= 0))) stop("wpop must be positive vector of length nrow(X)")
    if (any(w[predSampledAreas] != wpop[inSampleAreas])) stop("mismatch between w and wpop")
    if (!is.null(Narea)) {
      if (length(Narea) != M) stop("number of rows of Xpop does not match length of Narea")
      if (!is.null(names(Narea)) && !identical(names(Narea), rownames(Xpop))) warning("names of Narea not consistent with Xpop or sampled area names")
    }
    if (fpc) {
      if (is.null(Narea)) stop("area population sizes Narea required (unless fpc=FALSE)")
      if ((length(vpop) != M) || (any(vpop < 0))) stop("vpop must be a non-negative vector of length nrow(X)")
    }
  } else
    inSampleAreas <- predSampledAreas <- NULL

  # compute intermediate expressions to reduce computation time (of HB method)
  vx <- unname(as.matrix(crossprod(area.aggr, X)))
  vy <- as.vector(crossprod(area.aggr, y))
  if (equalv) {
    vnv <- vn / v
	  vxvb <- vx / vn
	  vyvb <- vy / vn
	  vv <- vn * v
	  sumlogv <- n * log(v)
	  XtX <- unname(as.matrix((crossprod(X) / v) + Omega0)) 
	  Xty <- as.vector((crossprod(X, y) / v) + Omega0 %*% beta0)
	  yyw <- as.numeric((sum(y^2) / v) + t(beta0) %*% Omega0 %*% beta0)
  } else {
    vnv <- as.vector(crossprod(area.aggr, 1/v))
    vxvb <- unname(as.matrix(crossprod(area.aggr, X/v)) / vnv)
    vyvb <- as.vector(crossprod(area.aggr, y/v)) / vnv
    vv <- as.vector(crossprod(area.aggr, v))
    sumlogv <- sum(log(v))
    XtX <- unname(as.matrix(crossprod(X/sqrt(v)) + Omega0))
    Xty <- as.vector(crossprod(X, y/v) + Omega0 %*% beta0)
    yyw <- as.numeric(sum(y^2 / v) + t(beta0) %*% Omega0 %*% beta0)
  }
  w.vnv <- w * vnv
  # number of degrees of freedom, last term is number of flat directions of beta prior
  dof <- nu0 + n - sum(eigen(Omega0, symmetric=TRUE, only.values=TRUE)$values == 0)
  nu0.s20 <- nu0 * s20
  if (M > 0) {
    fpc.factor <- 1L
    g1.0 <- rep.int(0, M)
    xdiff0 <- Xpop
    if (fpc) {  # requires area population sizes Narea
      g1.0 <- vpop/Narea
      g1.0[inSampleAreas] <- g1.0[inSampleAreas] - (vv[predSampledAreas]/Narea[inSampleAreas])/Narea[inSampleAreas]
      if (any(g1.0 < 0)) stop("vpop incompatible with v")
      fpc.factor <- 1 - (vn[predSampledAreas] / Narea[inSampleAreas])
	    if (any(fpc.factor < 0)) stop("area sample size exceeds area population size")
      xdiff0[inSampleAreas, ] <- xdiff0[inSampleAreas, ] - (vx[predSampledAreas, ] / Narea[inSampleAreas])
    }
    # NB fpc=FALSE can give smaller variances, because the O(1/N) term in g1 is ignored
  }

  lambdaREML <- lambdahat <- NULL

  if (!(method %in% c("synthetic", "survreg", "BLUP"))) {

    # log likelihood for lambda=sv2/se2
    LLH.lambda <- function(la) {
      out <- vector("double", length(la))
      out[la < 0] <- -Inf  # negative values impossible
      for (k in which(la >= 0)) {
        vnv.gamma <- la[k] * w.vnv; vnv.gamma <- vnv * vnv.gamma/(1 + vnv.gamma)
        H <- XtX - crossprod(vxvb * sqrt(vnv.gamma))
        xy <- Xty - colSums(vxvb * (vnv.gamma * vyvb))
        cholH <- chol(H); beta <- chol2inv(cholH) %*% xy  # instead of beta <- solve(H, xy)
        yKy <- yyw - sum(vnv.gamma * vyvb^2) - sum(xy * beta)
        out[k] <- compute.llh(la=la[k], logdetXX=2*sum(log(diag(cholH))), yKy=yKy)
      }
      out
    }

    # NB la must be scalar here
    compute.llh <- function(la, logdetXX, yKy) {
      if (nu0.s20 > yyw) {
        # for numerical stability; discard constant term dof * log(nu0.s20)
        temp <- dof * log(1 + (yKy/(nu0.s20)))
      } else {
        temp <- dof * log(nu0.s20 + yKy)
      }
      -0.5 * ( sum(log(1 + la * w.vnv)) + logdetXX + temp )
    }

    # Compute log likelihood at rough estimate of lambda to avoid over/underflow in computation of posterior probabilities for lambda
    if (is.null(lambda0)) {
      # pooled within area variance; discard areas with sample size 1
      if (equalv)
        se2.0 <- ((sum((vn - 1) * tapply(y, area, var)/v, na.rm=TRUE) / (n - m - sum(vn == 1))) + nu0.s20) / (1 + nu0)
      else
        se2.0 <- ((sum((vn - 1) * tapply(y, area, var)/tapply(v, area, mean), na.rm=TRUE) / (n - m - sum(vn == 1))) + nu0.s20) / (1 + nu0)
      sv2.0 <- max(var(vy/vn) - mean(v), var(vy/vn)/10)  # divide by 10 to allow for smaller residual variance ratio
      lambda0 <- sv2.0/se2.0
    }
    LLH.lambda0 <- LLH.lambda(lambda0)

    # REML estimate of lambda, do this anyway (i.e. also for HB or hybrid methods)
    global <- options(warn = -1)  # temporarily disable warnings
    opt <- NULL
    try(  # no need to stop if this fails
      opt <- optim(par=lambda0, LLH.lambda, method="BFGS", control=list(fnscale=-1, parscale=lambda0)), silent=TRUE
      #opt <- optim(par=lambda0, LLH.lambda, method="L-BFGS-B", lower=.Machine$double.eps^.5, upper=10*lambda0, control=list(fnscale=-1, parscale=lambda0))
    )
    options(global)
    if (!is.null(opt) && opt$convergence == 0)  # succesful optimization
      lambdaREML <- opt$par
    else {
      try({  # try another optimization method
        opt.tolerance <- min(.Machine$double.eps^0.25, lambda0/1e6)
        upperbound <- 1e4*lambda0
        opt <- optimize(LLH.lambda, c(0, upperbound), maximum=TRUE, tol=opt.tolerance)
        if (isTRUE(all.equal(opt$maximum, upperbound))) warning("reached upper bound of optimization interval")
        if (!silent && opt$maximum <= opt.tolerance) message("REML optimum below tolerance bound; may be zero")
        lambdaREML <- opt$maximum
      })
    }
    if ((method == "REML") && is.null(lambdaREML)) stop("REML optimization failed")
    if (!is.null(lambdaREML)) {
      if (lambdaREML <= .Machine$double.eps) lambdaREML <- 0
      if (!silent) message(sprintf("REML estimate of variance ratio: %.4g", lambdaREML))
      # change the starting 'value' for numerical integration (closer to REML; exactly REML turns out to slow down the integrations sometimes?)
      lambda0 <- 0.5 * (lambda0 + lambdaREML)
	}
    LLH.lambda0 <- LLH.lambda(lambda0)

    if (method != "REML") {

      # tables to hold lambda values and corr. (values prop. to) posterior probs.; stored for better performance
      lambda0.table <- NULL
      f0.table <- NULL

      # function to compute unnormalized posterior density for lambda
      # NOTE: rescale lambda (la) to focus on the most important integration region! (by factor lambda0)
      # The value of lambda0 (theoretically) just rescales the value of the integral by a factor 1/lambda0, yielding the same estimates
      # However, values of lambda0 far from the region of appreciable posterior mass may lead to numerical problems
      # (the numerical integrator may not be able to find that region in that case).
      post.lambda <- function(la) {
        post <- exp(LLH.lambda(la*lambda0) - LLH.lambda0)
        if (sum(is.infinite(post))) stop("Overflow; choose another starting value for lambda")
        post <- post * prior(la*lambda0)
        lambda0.table <<- c(lambda0.table, la)
        f0.table <<- c(f0.table, post)
        post
      }

      # function that returns x*f(x); used to compute posterior mean for lambda
      x.post.lambda <- function(la) {
        ind <- match(la, lambda0.table)
        if (any(is.na(ind)))
          f <- post.lambda(la)  # compute posterior if not yet computed for these x values
        else
          f <- f0.table[ind]  # otherwise refer to the table of already computed values
        la * lambda0 * f
      }

      # normalization constant
      if (!silent) cat("numerical integration of f(x) (normalization constant): ")
      numint0 <- integrate(post.lambda, lower=0, upper=Inf, rel.tol=.Machine$double.eps^0.5, ...)
      if (!silent) print(numint0, digits=max(3, getOption("digits") - 3))
      A <- numint0$value  # normalization constant
      if ((numint0$abs.error/A) > rel.int.tol)
        warning(sprintf("estimated relative numerical integration error %.4g for normalization constant", numint0$abs.error/A))

      # posterior mean for variance ratio lambda = s_v^2 / s_e^2
      if (!silent) cat("numerical integration of x*f(x): ")
      #numint1 <- integrate(x.post.lambda, lower=0, upper=Inf, abs.tol=A*lambda0/500, ...)
      numint1 <- integrate(x.post.lambda, lower=0, upper=Inf, rel.tol=.Machine$double.eps^0.5, ...)
      if (!silent) print(numint1, digits=max(3, getOption("digits") - 3))
      if ((numint1$abs.error/numint1$value) > rel.int.tol)
        warning(sprintf("estimated relative numerical integration error %.4g for posterior mean for between-area variance", numint1$abs.error/numint1$value))
      lambdahat <- numint1$value / A
      if (!silent) message(sprintf("posterior mean for variance ratio: %.4g", lambdahat))

      if (!silent) {
        posterior <- function(x) post.lambda(x/lambda0)/A
        curve(posterior, lambdahat/100, 5*lambdahat, main=expression(paste("posterior density for ", sigma[v]^2 / sigma[e]^2)), xlab=expression(sigma[v]^2 / sigma[e]^2), ylab="posterior density")
        # next two lines are for plotting the (scaled) prior for lambda in the same plot
        scaledprior <- function(x) prior(x) * 0.5 * posterior(lambdahat) / prior(lambdahat)
        curve(scaledprior, lambdahat/100, 5*lambdahat, lty=2, add=TRUE)
        abline(v=lambdahat, lty=3)  # add a vertical line at the posterior mean
        mtext("mean", side=1, at=lambdahat, cex=0.6, line=-0.5)
        if (!is.null(lambdaREML)) {
          abline(v=lambdaREML, lty=3, col=2)
          mtext("REML", side=1, at=lambdaREML, cex=0.6, col=2)
        }
      }
    }
  }

  # function to compute "BLUP" estimates
  # Note: this Bayesian BLUP already takes into account the variability of s_e^2!
  # Input:    la - variance ratio lambda
  # Output:   est - vector of BLUPs for area totals
  #	          mse - vector of MSEs of est
  #           synth - vector of synthetic components of est
  #           g1, g2 - components of mse (see e.g. Rao, 2003)
  #           gamma - vector of estimated variance ratios
  #	          beta - vector of estimated fixed effects
  #	          Vbeta - covariance matrix associated to beta
  #	          sigma2.hat - estimate of se2
  computeBLUP <- function(la) {

    gamma <- la * w.vnv; gamma <- gamma/(1 + gamma)

    # Although they usually do not affect the SAE estimates much, area effects are exluded from the fit
    # used to compute survey regression estimates.
    if (method == "survreg") {
      H <- XtX
      xy <- Xty
      yKy <- yyw
    } else {
      H <- XtX - crossprod(vxvb * sqrt(gamma * vnv))
      xy <- Xty - colSums(vxvb * (vyvb * gamma * vnv))
      yKy <- yyw - sum(gamma * vnv * vyvb^2)
    }

    cholH <- chol(H); P <- chol2inv(cholH)  # instead of P <- solve(H)
    beta <- P %*% xy
    yKy <- yKy - sum(xy * beta)
    sigma2.hat <- (nu0.s20 + yKy) / (dof - 2)  # BUP of s_e^2
    raneff <- gamma * (vyvb - as.numeric(vxvb %*% beta))

    # M=0 -> no prediction, only model fit
    if (M > 0) {
      synth <- as.numeric(Xpop %*% beta)
      if (fpc)
        synth[inSampleAreas] <- synth[inSampleAreas] + (vy - vx %*% beta)[predSampledAreas] / Narea[inSampleAreas]
      est <- synth
      est[inSampleAreas] <- est[inSampleAreas] + raneff[predSampledAreas] * fpc.factor
      # For general variance structure v, mBLUP cannot be decomposed as gamma*mSurvReg + (1-gamma)*mSynth,
      #   unless we ignore fpc; or simply define mSurvReg as (mBLUP - (1-gamma)*mSynth)/gamma

      g1 <- g1.0
      g1[-inSampleAreas] <- g1[-inSampleAreas] + wpop[-inSampleAreas] * la
      g1[inSampleAreas] <- g1[inSampleAreas] + (gamma/vnv)[predSampledAreas] * (fpc.factor^2)
      g1 <- sigma2.hat * g1

      xdiff <- xdiff0
      xdiff[inSampleAreas, ] <- xdiff[inSampleAreas, ] - fpc.factor * gamma[predSampledAreas] * vxvb[predSampledAreas, ]
      g2 <- sigma2.hat * rowSums(xdiff * (xdiff %*% P))  # rowSums(A * t(B)) is an efficient alternative for diag(A %*% B)	  
    }
    else
      synth <- est <- g1 <- g2 <- xdiff <- NULL

    list(est=est, mse=g1+g2, synth=synth, raneff=raneff, gamma=gamma, g1=g1, g2=g2,
         beta=beta, Vbeta=sigma2.hat*P, sigma2.hat=sigma2.hat, cholH=cholH, yKy=yKy, xdiff=xdiff)
  }

  # BLUP results for methods hybrid or REML; do this also for HB method (BLUP est and mse will be overwritten by HB estimates)
  lambda.plugin <- switch(method, synthetic=0, survreg=1e3*max(1/vnv), REML=lambdaREML, BLUP=lambda0, lambdahat)
  BLUP <- computeBLUP(lambda.plugin)

  names(BLUP$raneff) <- names(BLUP$gamma) <- sampledAreaNames
  BLUP$beta <- as.vector(BLUP$beta); names(BLUP$beta) <- colnames(X)
  dimnames(BLUP$Vbeta) <- list(colnames(X), colnames(X))
  if (M > 0)
    names(BLUP$est) <- names(BLUP$synth) <- rownames(Xpop)

  BLUP$Vraneff <- BLUP$gamma * (BLUP$sigma2.hat/vnv + rowSums(vxvb * (vxvb %*% BLUP$Vbeta)) * BLUP$gamma)
  BLUP$sigma2.se <- sqrt(2/(dof - 4)) * BLUP$sigma2.hat

  if ((M > 0) && full.cov) {
    BLUP$COV <- tcrossprod(BLUP$xdiff %*% t(chol(BLUP$Vbeta)))
    diag(BLUP$COV) <- diag(BLUP$COV) + BLUP$g1
  } else {
    BLUP$COV <- NULL
  }

  # compute model selection measures
  if (CV) {  # compute cross-validation measure; may take some time
    Xmin <- X - (BLUP$gamma * vxvb)[as.integer(area), ]  # Xmin generally not sparse
    BLUP$diagH <- ( (rowSums(Xmin * (Xmin %*% BLUP$Vbeta))/BLUP$sigma2.hat) + as.vector((BLUP$gamma/vnv)[as.integer(area)]) ) / v	
    BLUP$p.eff <- sum(BLUP$diagH)
    BLUP$res <- as.vector(X %*% BLUP$beta + BLUP$raneff[as.integer(area)])  # fitted values
    # compute fractions of fitted values exceeding the range of the data
    data.range <- range(y)
    BLUP$low.fitted <- sum(BLUP$res < data.range[1])/n
    BLUP$high.fitted <- sum(BLUP$res > data.range[2])/n
    BLUP$res <- y - BLUP$res  # residuals
    if (equalv)  # only compute R2 in the equal variance case
      BLUP$R2 <- 1 - (var(BLUP$res) / var(y))
    if (is.null(CVweights))
      CVweights <- 1 / (v + lambda.plugin)
    BLUP$CV <- mean( CVweights * (BLUP$res/(1 - BLUP$diagH))^2 ) / mean(CVweights)
    BLUP$llh.c <- -0.5 * ( n * log(2*pi*BLUP$sigma2.hat) + sumlogv + (1/BLUP$sigma2.hat) * sum((BLUP$res^2)/v) )
  }
  
  BLUP$cholH <- BLUP$yKy <- BLUP$xdiff <- NULL
  
  if ((M > 0) && (method == "HB")) {
    # initialize vectors for domain estimates
    estHB <- vector("double", M)
    mseHB <- vector("double", M)
    relErrM <- vector("double", M)
    relErrV <- vector("double", M)

    # tables to hold already computed quantities; for optimization purposes
    lambda.table <- NULL
    BLUP.table <- NULL
    mseBLUP.table <- NULL
    f.table <- NULL

    # function: compute.all() computes for all areas the required integrands for all variance ratios in la
    compute.all <- function(la) {

      BLUPtab <- matrix(NA_real_, length(la), M)
      mseBLUPtab <- matrix(NA_real_, length(la), M)
      f <- rep(NA_real_, length(la))

      for (k in 1:length(la)) {
        BLUP <- computeBLUP(la[k])

        BLUPtab[k,] <- BLUP$est
        mseBLUPtab[k, ] <- BLUP$mse  # diagonal of covariance matrix

        f[k] <- compute.llh(la=la[k], logdetXX=2*sum(log(diag(BLUP$cholH))), yKy=BLUP$yKy)
        f[k] <- prior(la[k]) * exp(f[k] - LLH.lambda0)

      }

      BLUP.table <<- rbind(BLUP.table, BLUPtab)
      mseBLUP.table <<- rbind(mseBLUP.table, mseBLUPtab)
      f.table <<- c(f.table, f)
    }

    # function: integr.ti
    # argument: la - variance ratio lambda (vectorized)
    #           i - domain number
    # returns:  BLUP_i(x) * f(x) at the points given in the vector x, where f(x) is  
    #           proportional to the posterior density for lambda and BLUP_i(x) is the BLUP estimate for domain i
    integr.ti <- function(la, i) {
      ind <- match(la, lambda.table)  # seems that either all or none of components of ind are NA
      if (any(is.na(ind))) {
        compute.all(la * lambda0)
        lambda.table <<- c(lambda.table, la)
        ind <- match(la, lambda.table)
      }
      BLUP.table[ind, i] * f.table[ind]
    }

    # function: integr.MSEi
    # argument: la - variance ratio lambda (vectorized)
    #           i - domain number
    #           post.mean - posterior mean for ith area total
    # returns: MSE_BLUP_i(x) * f(x) at the points given in the vector x, where f(x) is  
    #          proportional to the posterior density for lambda and MSE_BLUP_i(x) is the MSE of the BLUP estimate for domain i
    integr.MSEi <- function(la, i, post.mean) {
      ind <- match(la, lambda.table)
      if (any(is.na(ind))) {
        compute.all(la * lambda0)
        lambda.table <<- c(lambda.table, la)
        ind <- match(la, lambda.table)
      }
      (mseBLUP.table[ind, i] + (BLUP.table[ind, i] - post.mean)^2) * f.table[ind]
    }

    # compute HB estimates (posterior means and variances) for all areas
    abs.tol.M <- pmax(A*abs(BLUP$est)/200, .Machine$double.eps^0.5)
    abs.tol.V <- pmax(A*BLUP$mse/200, .Machine$double.eps^0.5)
    for (i in 1:M) {
      # compute posterior mean
      numint <- integrate(integr.ti, lower=0, upper=Inf, i=i, abs.tol=abs.tol.M[i], stop.on.error=FALSE, ...)
      relErrM[i] <- ifelse(abs(numint$value) > .Machine$double.eps, numint$abs.error/abs(numint$value), numint$abs.error)  # estimated integration errors
      if ((numint$message != "OK") || (relErrM[i] > rel.int.tol)) {  # another attempt with default tolerance values
        numint2 <- integrate(integr.ti, lower=0, upper=Inf, i=i, stop.on.error=FALSE, ...)
        if (numint2$abs.error < numint1$abs.error) {
          numint <- numint2
          relErrM[i] <- ifelse(abs(numint$value) > .Machine$double.eps, numint$abs.error/abs(numint$value), numint$abs.error)  # estimated integration errors
        }
      }
      if (numint$message != "OK") warning("error in numerical integration estimate ", i, ": ", numint$message)
      estHB[i] <- numint$value / A
      # compute posterior variance
      numint <- integrate(integr.MSEi, lower=0, upper=Inf, i=i, post.mean=estHB[i], abs.tol=abs.tol.V[i], stop.on.error=FALSE, ...)
      relErrV[i] <- ifelse(numint$value > .Machine$double.eps, numint$abs.error/numint$value, numint$abs.error)  # estimated integration errors
      if ((numint$message != "OK") || (relErrV[i] > rel.int.tol)) {  # another attempt with default tolerance values
        numint2 <- integrate(integr.MSEi, lower=0, upper=Inf, i=i, post.mean=estHB[i], stop.on.error=FALSE, ...)
        if (numint2$abs.error < numint1$abs.error) {
          numint <- numint2
          relErrV[i] <- ifelse(numint$value > .Machine$double.eps, numint$abs.error/numint$value, numint$abs.error)  # estimated integration errors
        }
      }
      if (numint$message != "OK") warning("error in numerical integration MSE ", i, ": ", numint$message)
      mseHB[i] <- numint$value / A
    }
    intM <- sum(relErrM > rel.int.tol, na.rm=TRUE)
    if (intM > 0) warning("number of inaccurate numerical integrations for posterior mean: ", intM)
    intV <- sum(relErrV > rel.int.tol, na.rm=TRUE)
    if (intV > 0) warning("number of inaccurate numerical integrations for posterior variance: ", intV)
    names(estHB) <- names(mseHB) <- names(relErrM) <- names(relErrV) <- rownames(Xpop)
  } else
    relErrV <- relErrM <- estHB <- mseHB <- NULL

  # output:
  if (method == "HB") {
    # (accounting for uncertainty in variance ratio would require M(M-1)/2 additional numerical integrations)
    if (!is.null(BLUP$COV)) {
	  f <- sqrt(mseHB/BLUP$mse)
	  BLUP$COV <- t(BLUP$COV * f) * f
	}
    BLUP$est <- estHB
    BLUP$mse <- mseHB
  }
  out <- c(BLUP, list(lambdaREML=lambdaREML, lambdahat=lambdahat, inSampleAreas=inSampleAreas, predSampledAreas=predSampledAreas,
           m=m, M=M, Narea=Narea, narea=vn, sampledAreaNames=sampledAreaNames, fpc=fpc,
           wpop=wpop,  # wpop needed to compute Vraneff for out-of-sample areas
           type="unit", method=method, relErrM=relErrM, relErrV=relErrV, aggregated=FALSE, benchmarked=FALSE))
  out$predAreaNames <- rownames(Xpop)
  names(out$narea) <- sampledAreaNames
  if (CV || keep.data)
    out <- c(out, list(y=y))
  if (keep.data) {
    out$xs <- vx / vn
    dimnames(out$xs) <- list(sampledAreaNames, colnames(Xpop))
    out <- c(out, list(X=X, Xp=Xpop, w=w, v=v, vpop=vpop, area=area,
             prior=prior, beta0=beta0, Omega0=Omega0, nu0=nu0, s20=s20))
  }
  class(out) <- "sae"
  out
}


#' Compute area-level cross-validation measure for sae objects.
#'
#' This function computes a cross-validation measure defined at the area level.
#' It can be used, for example, to compare the predictive ability of area and unit-level models.
#' The code is based in part on that of \code{cv.glm} from package \pkg{boot}.
#'
#' @export
#' @param sae object of class sae, resulting from a call to \code{fSAE}, \code{fSAE.Area}, or \code{fSAE.Unit}.
#' @param weight if \code{TRUE}, use weights inversely proportional to the MSEs of y - yhat in the cost function.
#' @param cost cost function to be used. Defaults to a quadratic cost function.
#' @param K K in K-fold cross-validation. Specifies in how many parts the dataset should be divided.
#' @param method method used to refit the model. One of "HB", "hybrid" (default) or "REML", in the order of slow to fast.
#' @param seed random seed used in selecting groups of areas to leave out in K-fold cross-validation.
#' @return The computed area-level cross-validation measure.
#' @example examples/CVarea.R
CVarea <- function(sae, weight=TRUE, cost = function(y, yhat, w) sum(w * (y - yhat)^2)/sum(w), K=10, method="hybrid", seed) {
  if (!is(sae, "sae") || sae$type == "direct") stop("not a model-based sae object")
  if (is.null(sae[["X"]])) stop("sae object does not contain input data; regenerate it using keep.data=TRUE")

  n <- length(sae$sampledAreaNames)  # number of areas in the sample

  if ((K > n) || (K <= 1)) stop("K outside allowable range")
  K0 <- K
  K <- round(K)
  kvals <- unique(round(n/(1L:floor(n/2))))
  temp <- abs(kvals - K)
  if (!any(temp == 0)) 
    K <- kvals[temp == min(temp)][1L]
  if (K != K0)
    warning("K has been set to ", K)
  f <- ceiling(n/K)  # max number of areas to leave out
  if (!missing(seed))
    set.seed(seed)
  s <- sample(rep.int(1L:K, f), n)
  ms <- max(s)

  CVp <- vector("double", n)  # predictions
  CVr <- vector("double", n)  # reference values
  CVw <- vector("double", n)  # weights

  # CV to measure deviation from survey regression estimates
  if (sae$type == "area") {  # area-level model
    sreg <- list(est=sae$y, mse=sae$v)  # survey regression estimates already computed for area level
    dimnames(sae$X) <- dimnames(sae$Xp[sae$sampledAreaNames, ])
  } else  # unit-level model
    sreg <- fSurvReg(sae$y, sae$X, sae$area, sae$Narea, sae$Xp, v=sae$v, vpop=sae$vpop)
  index <- 1
  for (i in seq_len(ms)) {
    j.in <- sae$sampledAreaNames[s != i]
    j.out <- sae$sampledAreaNames[s == i]
    m <- length(j.out)
	if (m > 0) {
	  if (sae$type == "area") {
        r <- fSAE.Area(sreg$est[j.in], sreg$mse[j.in], sae$Xp, sae$X, Narea=sae$Narea, method=method, silent=TRUE, full.cov=FALSE, CV=FALSE,
                       w=sae$w[sae$area %in% j.in], wpop=sae$wpop, beta0=sae$beta0, Omega0=sae$Omega0, nu0=sae$nu0, s20=sae$s20, prior=sae$prior)
      } else {  # unit-level model
        r <- fSAE.Unit(sae$y[sae$area %in% j.in], sae$X[sae$area %in% j.in, ], as.factor(sae$area[sae$area %in% j.in]), sae$Narea, sae$Xp, method=method, silent=TRUE, full.cov=FALSE, CV=FALSE,
                       fpc=sae$fpc, v=if (length(sae$v) == 1) sae$v else sae$v[sae$area %in% j.in], vpop=sae$vpop, w=sae$w[sae$sampledAreaNames %in% j.in], wpop=sae$wpop, beta0=sae$beta0, Omega0=sae$Omega0, nu0=sae$nu0, s20=sae$s20, prior=sae$prior)
      }
      CVp[index:(index + m - 1)] <- r$est[j.out]
      CVr[index:(index + m - 1)] <- sreg$est[j.out]
      CVw[index:(index + m - 1)] <- ifelse(weight, 1/(sreg$mse[j.out] + r$mse[j.out]), rep.int(1L, m))
      index <- index + m
	}
  }
  cost(CVr, CVp, CVw)
}


#' Compute aggregates of small area estimates and MSEs.
#'
#' @export
#' @param x sae object.
#' @param R aggregation matrix, r x M matrix where M is the number of areas and r the number of aggregate areas; default is aggregation over all areas.
#' @return Object of class \code{sae} with aggregated small area estimates and MSEs.
#' @seealso \code{\link{sae-class}}
#' @example examples/aggr.R
# NB: gamma, raneff, Vraneff, diagH, ... are associated with the model fit and are unaffected by aggr()
aggr <- function(x, R) {
  if (is.null(x$Narea)) stop("x requires Narea component containing area population sizes")
  if (missing(R)) R <- Matrix(1, nrow=x$M)
  if (!is(R, "Matrix")) R <- Matrix(R)
  if (nrow(R) != x$M) stop("matrix R has incorrect dimensions")
  Naggr <- as.vector(crossprod(R, x$Narea))
  names(Naggr) <- dimnames(R)[[2]]  # use column names of Naggr, if any
  RN <- x$Narea * R
  x$est <- as.vector(crossprod(RN, x$est)) / Naggr
  if (!is.null(x$COV)) {
    x$COV <- as.matrix(crossprod(RN, x$COV) %*% RN) / Naggr^2
    x$mse <- diag(x$COV)
  } else {
    if (x$type != "direct") warning("no covariance matrix available; zero covariances assumed")
    x$mse <- as.vector(crossprod(RN, x$Narea * x$mse)) / Naggr^2
  }
  if (!is.null(x$g1)) x$g1 <- as.vector(crossprod(RN, x$Narea * x$g1)) / Naggr^2
  if (!is.null(x$g2)) x$g2 <- x$mse - x$g1

  x$aggregated <- TRUE
  x
}

#' Benchmark small area estimates.
#'
#' Benchmark small area estimates to agree with given totals at aggregate levels.
#'
#' @export
#' @param x sae object to be benchmarked.
#' @param R restriction matrix, r x M matrix where r is the number of restrictions and M the number of areas; default is a single constraint on the population total.
#'    Note that R acts on the vector of area population totals, not the vector of means.
#' @param rhs r-vector of benchmark totals corresponding to the restrictions in the rows of R.
#' @param mseMethod if "no", MSEs are not updated, if "exact", constraints are treated as exact identities, and
#'    if "model", the squared differences between original and benchmarked estimates are added to the MSEs.
#  mseMethod="model" is like assuming that the model is true (see You, Rao, Dick paper)
#' @return An object of class \code{sae} with adjusted estimates.
#' @seealso \code{\link{sae-class}}
#' @references
#'   Y. You, J.N.K. Rao and P. Dick (2004). Benchmarking Hierarchical Bayes Small Area Estimators
#'     in the Canadian Census Undercoverage Estimation. Statistics in Transition 6(5), 631-640.
#' @example examples/bench.R
bench <- function(x, R, rhs, mseMethod="no") {
  if (is.null(x$Narea)) stop("x requires Narea component containing area population sizes")
  if (missing(R)) R <- Matrix(1, nrow=x$M)
  if (!is(R, "Matrix")) R <- Matrix(R)
  if (nrow(R) != x$M) stop("matrix R has incorrect dimensions")
  if (ncol(R) != length(rhs)) stop("matrix R and rhs incompatible")
  if (is.null(x$COV))
    V <- (x$Narea^2) * x$mse
  else
    V <- symmpart(t(x$COV * x$Narea) * x$Narea)
  res <- P.update(x$Narea * x$est, V, R, rhs, 0, cov = (mseMethod=="exact"))
  if (mseMethod != "no") {
    if (mseMethod == "exact") {
      x$mse <- diag(res$V) / x$Narea^2
      if (!is.null(x$COV))
        x$COV <- symmpart(t(res$V / x$Narea) / x$Narea)
    } else {  # add squared differences as a squared bias term
      mse0 <- x$mse
      x$mse <- x$mse + (x$est - (res$M/x$Narea))^2
      if (!is.null(x$COV)) {
        f <- x$mse/mse0
        x$COV <- t(x$COV * f) * f
      }
    }
  }
  x$est <- res$M / x$Narea
  x$benchmarked <- TRUE
  x
}

#' Compute unit weights underlying the small area estimates or their aggregate.
#'
#' The small area estimates can be interpreted as weighted sums of the response variable. This function computes the weights
#' corresponding to the aggregated small area estimates or the weights corresponding to a specific small area estimate. The
#' weights applied to the response variable need not exactly reproduce the Hierarchical Bayes estimate since the latter
#' is averaged over the posterior distribution for the variance ratio whereas the weights are evaluated at the posterior mean.
#' Under the default prior for the fixed effects, the weights applied to the design matrix reproduce the corresponding population numbers.
#'
#' @export
#' @param x sae object.
#' @param areaID if left unspecified (\code{NULL}), weights corresponding to the overall (aggregated) estimate are returned.
#'   Otherwise weights that reproduce the estimate for a specific area are returned.
#' @param forTotal if \code{FALSE} weights will be divided by the population size.
#' @return An object of class \code{weights}.
#' @seealso \code{\link{summary.weights}}, \code{\link{plot.weights}}
#' @example examples/uweights.R
uweights <- function(x, areaID=NULL, forTotal=FALSE) {
  if (is.null(x[["X"]])) stop("sae object does not contain input data; regenerate it using keep.data=TRUE")
  if (any(x$beta0 != 0)) stop("zero prior mean beta0 required")
  if (x$aggregated || x$benchmarked) stop("not available for benchmarked/aggregated sae objects")
  areavar <- as.character(x$area)
  if (is.null(areaID))
    w <- as.vector(1 - x$gamma[areavar] + x$gamma[areavar]*(x$Narea[areavar]/x$narea[areavar]))
  else {
    areaID <- as.character(areaID)
	if (!(areaID %in% names(x$Narea))) stop("unrecognized areaID")
    w <- ifelse(x$area == areaID, 1, 0) * (1 - x$gamma[areavar] + x$gamma[areavar]*(x$Narea[areavar]/x$narea[areavar]))
  }
  Xp.adj <- x$Narea * x$Xp
  Xp.adj[x$sampledAreaNames, ] <- Xp.adj[x$sampledAreaNames, ] - x$narea * x$xs - x$gamma * (x$Narea[x$sampledAreaNames] - x$narea) * x$xs
  if (is.null(areaID))
    vec <- colSums(Xp.adj)
  else
    vec <- Xp.adj[areaID, ]
  temp <- as.vector(x$X %*% crossprod(x$Vbeta, vec / x$sigma2.hat)) / x$v
  area.aggr <- aggrMatrix(x$area)
  if (length(x$v) == 1)
    vnv <- x$narea
  else
    vnv <- as.vector(crossprod(area.aggr, 1/x$v))
  vxvb <- unname(as.matrix(crossprod(area.aggr, x$X/x$v)) / vnv)
  vtempvb <- as.vector(crossprod(area.aggr, temp)) / vnv
  w <- w + temp - as.vector((vtempvb * x$gamma)[areavar]) / x$v
  if (!forTotal) {
    if (is.null(areaID))
      w <- w / sum(x$Narea)
    else
      w <- w / x$Narea[areaID]
  }
  class(w) <- "weights"
  w
}


###################
# Output functions

#' Print method for objects of class sae.
#'
#' @S3method print sae
#' @method print sae
#' @param x object of class \code{sae}.
#' @param digits number of digits to display.
#' @param ... additional arguments passed to \code{print.default}.
print.sae <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  varcomp <- c("posterior mean sv2" = x$sv2hat, "REML estimate sv2" = x$sv2REML, "posterior mean sv2/se2" = x$lambdahat, "REML estimate sv2/se2" = x$lambdaREML)
  if (!is.null(x$call)) cat("Call: ", deparse(x$call), "\n\n", sep = "")
  cat("Summaries:\n")
  summ <- rbind("gamma weights" = summary(wDirect(x)), "random effects" = summary(raneff(x)))
  aggrs <- NULL
  if (!is.null(x$est)) {
    summ <- rbind(estimates = summary(x$est), "standard errors" = summary(sqrt(x$mse)), "rel. std. err." = summary(relSE(x)), summ)
    if (x$M)
      aggrs <- list(mean = sum(x$Narea * x$est)/sum(x$Narea), total = sum(x$Narea * x$est))
  }
  print.default(format(summ, digits = digits), print.gap = 2, quote = FALSE, ...)
  if (x$benchmarked) cat("\nEstimates have been benchmarked.\n")
  if (x$aggregated) cat("\nEstimates have been aggregated.\n")
  if (!is.null(aggrs)) {
    cat("\nPopulation aggregates:\n")
    print.default(format(aggrs, digits = digits), print.gap = 2, quote = FALSE, ...)
  }
  varcomp <- c(varcomp, "se2 estimate" = se2(x))
  if (length(varcomp)) {
    cat("\n")
    print.default(format(varcomp, digits = digits), print.gap = 2, quote = FALSE, ...)
  }
  if (length(coef(x))) {
    cat("\nCoefficients:\n")
    coefs <- coef(x)
    if (length(vcov(x))) {
      coefs <- cbind(coefs, sqrt(diag(vcov(x))), coefs/sqrt(diag(vcov(x))))
      dimnames(coefs)[[2]] <- c("coefficient", "standard error", "t value")
    }
    printCoefmat(coefs, digits = digits - 1)
  }
  modsel <- list(CV = x$CV, "p.eff" = x$p.eff, cAIC = cAIC(x))
  if (length(x$CV)) {
    cat("\nModel selection measures:\n")
    print.default(format(modsel, digits = digits), print.gap = 2, quote = FALSE, ...)
    if (length(x$R2))
      message("R-squared: ", round(x$R2, 5))
    if (!is.null(x$low.fitted) && x$low.fitted > 0)
      message("Fraction of fitted values below lowest data value: ", round(x$low.fitted, 5))
    if (!is.null(x$high.fitted) && x$high.fitted > 0)
      message("Fraction of fitted values above highest data value: ", round(x$high.fitted, 5))
  }
}

#' Plot method for objects of class sae.
#'
#' This function calls \code{coefplot} from package \pkg{arm} to plot small area estimates and standard errors.
#' Multiple sets of estimates can be compared using this plot. The default ordering of the estimates is by their
#' area population sizes.
#'
#' @S3method plot sae
#' @method plot sae
#' @param ... list of \code{sae} objects (or lists with at least a component "est" (and optionally "mse")).
#'   Other components not matching other arguments are passed to the lower-level plotting functions.
#' @param sort.by vector by which to sort the vectors to be plotted; defaults to population size vector of the first sae object.
#' @param decreasing if \code{TRUE}, sort in decreasing order (default).
#' @param index vector of indices of the selected areas to be plotted.
#' @param maxrows maximum number of rows in a column.
#' @param maxcols maximum number of columns of estimates on a page.
#' @param labels optional vector of lables; default is the vector of area names of the first sae object.
#' @param type "sae" for small area estimates (default), "coef" for coefficients, "raneff" for random effects.
#' @param ylab label for y-axis.
#' @param main main title.
#' @param offset space used between plots of multiple estimates for the same area.
#' @param cex.var the fontsize of the variable names, default=0.9.
#' @param mar a numerical vector of the form c(bottom, left, top, right) which gives the number of lines of margin to be specified on the four sides of the plot. The default is c(0,1,5,0).
#' @param CI confidence interval, default is 1, which will plot plus and minus 1 standard errors. If CI=2, plots plus and minus 2 standard errors instead.
#' @importFrom arm coefplot coefplot.default
plot.sae <- function(..., sort.by=NULL, decreasing=FALSE, index=NULL,
              maxrows=50, maxcols=6, labels=NULL, type="sae",
			  ylab = if (type == "coef") "coef" else "area",
			  main = switch(type, coef="coefficients", raneff="random effects", "SAE estimates"),
              offset=0.2, cex.var=0.9, mar=c(0,1,5,0), CI=1) {

  dotargs <- list(...)
  toplot <- sapply(dotargs, function(obj) is.list(obj) && !is.null(EST(obj, type)))
  x <- dotargs[toplot]
  if (!length(x)) stop("nothing to plot")
  est.plot <- lapply(x, EST, type=type)
  se.plot <- lapply(x, SE, type=type)
  se.plot <- lapply(se.plot, function(obj) {if (is.null(obj)) obj <- 0 * est.plot[[1]]; obj})
  grpar <- dotargs[!toplot]  # other graphical parameters
  if (!length(grpar)) grpar <- NULL

  if (is.null(sort.by))
    if (!(type %in% c("coef", "raneff")) && !is.null(x[[1]]$Narea) && (length(x[[1]]$Narea) == length(est.plot[[1]])))
      sort.by <- 1/x[[1]]$Narea  # default: order from large to small areas 
    else
      sort.by <- 1:length(est.plot[[1]])

  if (is.null(index)) {
    M <- length(est.plot[[1]])
    o <- order(sort.by, decreasing=decreasing)
  } else {
    M <- length(index)
    o <- index[order(sort.by[index], decreasing=decreasing)]
  }

  if (!is.null(labels))
    lab <- labels[o]
  else {
    if (!is.null(names(est.plot[[1]])))
	  lab <- names(est.plot[[1]][o])
	else
	  lab <- switch(type, coef=1:length(est.plot[[1]]), raneff=x$sampledAreaNames, names(x[[1]]$Narea))[o]
  }
  if (is.null(lab)) lab <- seq_along(o)

  maxrows <- min(maxrows, M)

  cols <- ceiling(M/maxrows)
  pages <- ceiling(cols/maxcols)

  for (page in 1:pages) {
    colrange <- ((page - 1) * maxcols + 1):min((page * maxcols), cols)
    dev.new()
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mfrow = c(1, min(maxcols, cols)))
    for (co in colrange) {
      par(mfg=c(1, co - (colrange[1] - 1)))  # force plot in the intended column
      rowrange <- ((co - 1) * maxrows + 1):(co * maxrows)
      offsets <- offset * (seq_along(x) - 1)  # - ((length(x) + 1) %/% 2))
      for (i in seq_along(x)) {
        if (i==1) {
          do.call(coefplot, c(list(est.plot[[i]][o[rowrange]], se.plot[[i]][o[rowrange]], varnames=lab[rowrange],
            cex.var=cex.var, main=main, ylab=ylab, offset=offsets[i], mar=mar, CI=CI), grpar))
        } else {
          coefplot(est.plot[[i]][o[rowrange]], se.plot[[i]][o[rowrange]], col.pts=i, add=TRUE, offset=offsets[i], mar=mar, CI=CI)
        }
      }
    }
  }

}

#' Summary method for objects of class \code{weights}.
#'
#' @S3method summary weights
#' @method summary weights
#' @param object object of class \code{weights} as returned by function \code{\link{uweights}}.
#' @param ... not used.
#' @seealso \code{\link{uweights}}, \code{\link{plot.weights}}
summary.weights <- function(object, ...) {
  deft <- sqrt(deff.weights(object))
  object <- unclass(object)
  cat("\n")
  print(summary(object))
  cat("\n")
  cat(sprintf("                       number:  %i \n", length(object) ))
  cat(sprintf("                          sum:  %.1f \n", sum(object) ))
  cat(sprintf("   number of negative weights:  %i \n", sum(object < 0) ))
  cat(sprintf("          number of 0 weights:  %i \n", sum(object == 0) ))
  cat(sprintf("number of weights > 0 and < 1:  %i \n", sum(object > 0 & object < 1) ))
  cat(sprintf("          number of weights 1:  %i \n", sum(object == 1) ))
  cat(sprintf("          sqrt(design-effect):  %.4f \n\n", deft ))
}

#' Plot method for objects of class \code{weights}.
#'
#' @S3method plot weights
#' @method plot weights
#' @param x object of class \code{weights} as returned by function \code{\link{uweights}}.
#' @param log whether to log-tranform the weights.
#' @param main main title of plot.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param col colour.
#' @param ... additional arguments passed to \code{\link{hist}}.
#' @seealso \code{\link{uweights}}, \code{\link{summary.weights}}
plot.weights <- function(x, log=FALSE, main="Distribution of weights",
                         xlab=if (log) "log(weight)" else "weight", ylab="frequency", col="cyan", ...) {
  x <- unclass(x)
  if (log) x <- log(x)
  h <- hist(x, breaks="Scott", main=main, xlab=xlab, ylab=ylab, col=col, ...)
  boxplot(x, horizontal=TRUE, add=TRUE, at=-max(h$counts)/50, axes=TRUE)
}


######################
# Extractor functions

#' S3 class for the fitted model and SAE outcomes.
#'
#' Functions \code{\link{fSAE}}, \code{\link{fSurvReg}}, \code{\link{fSAE.Area}} and \code{\link{fSAE.Unit}}
#' return an object of class \code{sae}. It contains information on the model fit as well as the
#' small area estimates, error estimates and a few model selection measures.
#' The functions listed below extract the main components from an object of class \code{sae}.
#' \describe{
#'  \item{\code{EST(x, type="sae")}}{return the vector of small area estimates of \code{sae} object x. Alternatively,
#'   with \code{type} "coef" or "raneff" fixed or random effect estimates are returned.}
#'  \item{\code{MSE(x, type="sae")}}{return the vector of mean squared errors of \code{sae} object x. Alternatively,
#'   with \code{type} "coef" or "raneff" MSEs of fixed or random effects are returned.}
#'  \item{\code{SE(x, type="sae")}}{extract standard errors, i.e. square roots of the MSEs.}
#'  \item{\code{relSE(x, type="sae")}}{extract relative standard errors.}
#'  \item{\code{COV(x)}}{extract the covariance matrix for the small area estimates.}
#'  \item{\code{COR(x)}}{extract the correlation matrix for the small area estimates.}
#'  \item{\code{coef(x)}}{\code{coef} method for \code{sae} objects; returns vector of fixed effects.}
#'  \item{\code{vcov(x)}}{\code{vcov} method for \code{sae} objects; returns covariance matrix for fixed effects.}
#'  \item{\code{raneff(x, pop)}}{return vector of random effects. If \code{pop=TRUE} returns a vector for predicted areas (zero for out-of-sample areas), otherwise a vector for sampled areas.}
#'  \item{\code{raneff.se(x, pop)}}{return vector of standard errors for random effects.}
#'  \item{\code{residuals(x)}}{\code{residuals} method for \code{sae} objects; returns a vector of residuals.}
#'  \item{\code{fitted(x)}}{\code{fitted} method for \code{sae} objects; returns a vector of fitted values.}
#'  \item{\code{se2(x)}}{extracts within-area variance estimate.}
#'  \item{\code{sv2(x)}}{extracts between-area variance estimate.}
#'  \item{\code{wDirect(x, pop)}}{extract vector of weights of the survey regression components in the small area estimates. If \code{pop=TRUE} returns a vector for predicted areas (zero for out-of-sample areas), otherwise a vector for sampled areas.}
#'  \item{\code{synthetic(x)}}{extract vector of synthetic estimates.}
#'  \item{\code{CV(x)}}{extract leave-one-out cross-validation measure.}
#'  \item{\code{cAIC(x)}}{extract conditional AIC measure.}
#'  \item{\code{R2(x)}}{extract unit-level R-squared goodness-of-fit measure.}
#' }
#' Other components include
#' \describe{
#'  \item{\code{relErrM,relErrV}}{relative numerical integration errors in estimates and MSEs, for \code{method} "HB".}
#' }
#'
#' @name sae-class
#' @aliases EST MSE SE COV COR relSE coef.sae vcov.sae raneff raneff.se residuals.sae fitted.sae se2 sv2 wDirect synthetic CV cAIC R2
#' @example examples/sae-class.R
NULL

#' @export
#' @noRd
EST <- function(x, type="sae") {
  switch(type,
    coef = coef(x),
    raneff = raneff(x),
    x$est
  )
}

#' @export
#' @noRd
MSE <- function(x, type="sae") {
  switch(type,
    coef = diag(x$Vbeta),
    raneff = x$Vraneff,
    x$mse
  )
}

#' @export
#' @noRd
SE <- function(x, type="sae") {
  out <- MSE(x, type)
  if (is.null(out))
    NULL
  else
    sqrt(out)
}

#' @export
#' @noRd
COV <- function(x) {
  x$COV
}

#' @export
#' @noRd
COR <- function(x) {
  out <- COV(x)
  if (!is.null(out)) {
    f <- 1/sqrt(x$mse)
    out <- t(out * f) * f
  }
  out
}

#' @export
#' @noRd
relSE <- function(x, type="sae") {
  SE(x, type)/EST(x, type)
}

#' @S3method coef sae 
#' @noRd
coef.sae <- function(object, ...) {
  object$beta
}

#' @S3method vcov sae
#' @noRd
vcov.sae <- function(object, ...) {
  object$Vbeta
}

#' @export
#' @noRd
raneff <- function(x, pop=(x$M > 0)) {
  if (pop) {  # predicted areas; zero for out-of-sample areas
    out <- rep(0, x$M)
    names(out) <- x$predAreaNames
    out[x$inSampleAreas] <- x$raneff[x$predSampledAreas]
  } else  # sampled areas
    out <- x$raneff
  out
}

#' @export
#' @noRd
raneff.se <- function(x, pop=(x$M > 0)) {
  if (pop) {  # predicted areas
    out <- x$wpop * sv2(x)
    names(out) <- x$predAreaNames
    out[x$inSampleAreas] <- x$Vraneff[x$predSampledAreas]
  } else  # sampled areas
    out <- sqrt(x$Vraneff)
  sqrt(out)
}

#' @S3method residuals sae
#' @noRd
residuals.sae <- function(object, ...) {
  if (is.null(object$res)) warning("use CV=TRUE to compute residuals")
  object$res
}

#' @S3method fitted sae
#' @noRd
fitted.sae <- function(object, ...) {
  object$y - residuals(object)
}

#' @export
#' @noRd
se2 <- function(x) {
  if (class(x) != "sae") stop("require sae object as input")
  x$sigma2.hat
}

#' @export
#' @noRd
sv2 <- function(x) {
  if (class(x) != "sae") stop("require sae object as input")
  switch(x$method,
    survreg = Inf,
    synthetic = 0,
    REML = x$lambdaREML * se2(x),
    x$lambdahat * se2(x)
  )
}

#' @export
#' @noRd
wDirect <- function(x, pop=(x$M > 0)) {
  if (pop) {  # predicted areas; zero for out-of-sample areas
    out <- rep(0, x$M)
    names(out) <- x$predAreaNames
    out[x$inSampleAreas] <- x$gamma[x$predSampledAreas]
  } else  # sampled areas
    out <- x$gamma
  out
}

#' @export
#' @noRd
synthetic <- function(x) {
  switch(x$method,
    survreg = rep(0, x$M),
    x$synth
  )
}

#' @export
#' @noRd
CV <- function(x) {
  x$CV
}

#' @export
#' @noRd
cAIC <- function(x) {
  - 2 * x$llh.c + 2 * x$p.eff
}

#' @export
#' @noRd
R2 <- function(x) {
  x$R2
}


##################
# Other functions

#' alphabetically order names in interactions
#' @noRd
orderInteractions <- function(names) {
  sapply( strsplit(names, ":"), function(s) {paste(sort(s), collapse=":")} )
}

#' compute design effect due to (unequal) weights
#' @noRd
deff.weights <- function(w) {
  1 + (var(w)/(mean(w))^2)
}

#' Utility function to construct a sparse aggregation matrix from a factor.
#'
#' @noRd
#' @param fac factor variable.
#' @param facnames whether the factor levels should be used as column names for the aggregation matrix.
#' @return Sparse aggregation matrix of class \code{\link[=dgCMatrix-class]{dgCMatrix}}.
aggrMatrix <- function(fac, facnames=FALSE) {
  fac <- as.factor(fac)
  if ((length(fac) == nlevels(fac)) && all(fac == levels(fac)))
    Maggr <- Diagonal(length(fac))
  else
    Maggr <- sparseMatrix(i=1:length(fac), j=as.integer(fac), x=rep.int(1L, length(fac)), dims=c(length(fac), nlevels(fac)))
  if (facnames)
    dimnames(Maggr) <- list(NULL, levels(fac))
  Maggr
}

#' Bayes update of normal prior to posterior under normal data on linear combinations
#' @noRd
P.update <- function(priorM, priorV, Dmat, d, Omega=0, cov=TRUE) {
  # if priorV is a vector assume a diagonal covariance matrix
  if (is.vector(priorV)) {
    if (length(priorV) != length(priorM))
      stop("lengths of priorM and priorV incompatible")
    priorV <- Diagonal(x=priorV)
  }
  DV <- crossprod(Dmat, priorV)
  DVD <- DV %*% Dmat
  DVD <- symmpart(DVD + Omega)
  invDVD <- chol2inv(chol(DVD))
  G <- crossprod(DV, invDVD)
  M <- as.vector(priorM + G %*% (d - crossprod(Dmat, priorM)))
  if (cov)
    V <- priorV - G %*% DV
  else
    V <- NULL
  list(M=M, V=V)
}
