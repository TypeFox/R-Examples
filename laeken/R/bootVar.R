# ---------------------------------------
# Author: Andreas Alfons
#         Vienna University of Technology
# ---------------------------------------

## TODO: support estimators based on semiparametric outlier detection
## FIXME: do not use 'p' as argument name for function passed to 'boot'

#' Bootstrap variance and confidence intervals of indicators on social exclusion
#' and poverty
#'
#' Compute variance and confidence interval estimates of indicators on social
#' exclusion and poverty based on bootstrap resampling.
#'
#' @param inc either a numeric vector giving the equivalized disposable income,
#' or (if \code{data} is not \code{NULL}) a character string, an integer or a
#' logical vector specifying the corresponding column of \code{data}.
#' @param weights optional; either a numeric vector giving the personal sample
#' weights, or (if \code{data} is not \code{NULL}) a character string, an
#' integer or a logical vector specifying the corresponding column of
#' \code{data}.
#' @param years optional; either a numeric vector giving the different years of
#' the survey, or (if \code{data} is not \code{NULL}) a character string, an
#' integer or a logical vector specifying the corresponding column of
#' \code{data}.  If supplied, values are computed for each year.
#' @param breakdown optional; either a numeric vector giving different domains,
#' or (if \code{data} is not \code{NULL}) a character string, an integer or a
#' logical vector specifying the corresponding column of \code{data}.  If
#' supplied, the values for each domain are computed in addition to the overall
#' value.
#' @param design optional; either an integer vector or factor giving different
#' strata for stratified sampling designs, or (if \code{data} is not
#' \code{NULL}) a character string, an integer or a logical vector specifying
#' the corresponding column of \code{data}.  If supplied, this is used as
#' \code{strata} argument in the call to \code{\link[boot]{boot}}.
#' @param cluster optional; either an integer vector or factor giving different
#' clusters for cluster sampling designs, or (if \code{data} is not
#' \code{NULL}) a character string, an integer or a logical vector specifying
#' the corresponding column of \code{data}.
#' @param data an optional \code{data.frame}.
#' @param indicator an object inheriting from the class \code{"indicator"} that
#' contains the point estimates of the indicator (see \code{\link{arpr}},
#' \code{\link{qsr}}, \code{\link{rmpg}} or \code{\link{gini}}).
#' @param R a numeric value giving the number of bootstrap replicates.
#' @param bootType a character string specifying the type of bootstap to be
#' performed.  Possible values are \code{"calibrate"} (for calibration of the
#' sample weights of the resampled observations in every iteration) and
#' \code{"naive"} (for a naive bootstrap without calibration of the sample
#' weights).
#' @param X if \code{bootType} is \code{"calibrate"}, a matrix of calibration
#' variables.
#' @param totals numeric; if \code{bootType} is \code{"calibrate"}, this gives
#' the population totals.  If \code{years} is \code{NULL}, a vector should be
#' supplied, otherwise a matrix in which each row contains the population totals
#' of the respective year.  If this is \code{NULL} (the default), the population
#' totals are computed from the sample weights using the Horvitz-Thompson
#' estimator.
#' @param ciType a character string specifying the type of confidence
#' interval(s) to be computed.  Possible values are \code{"perc"}, \code{"norm"}
#' and \code{"basic"} (see \code{\link[boot]{boot.ci}}).
#' @param alpha a numeric value giving the significance level to be used for
#' computing the confidence interval(s) (i.e., the confidence level is \eqn{1 -
#' }\code{alpha}), or \code{NULL}.
#' @param seed optional; an integer value to be used as the seed of the random
#' number generator, or an integer vector containing the state of the random
#' number generator to be restored.
#' @param na.rm a logical indicating whether missing values should be removed.
#' @param gender either a numeric vector giving the gender, or (if \code{data}
#' is not \code{NULL}) a character string, an integer or a logical vector
#' specifying the corresponding column of \code{data}.
#' @param method a character string specifying the method to be used (only for
#' \code{\link{gpg}}).  Possible values are \code{"mean"} for the mean, and
#' \code{"median"} for the median.  If weights are provided, the weighted mean
#' or weighted median is estimated.
#' @param \dots if \code{bootType} is \code{"calibrate"}, additional arguments
#' to be passed to \code{\link{calibWeights}}.
#'
#' @return An object of the same class as \code{indicator} is returned.  See
#' \code{\link{arpr}}, \code{\link{qsr}}, \code{\link{rmpg}} or
#' \code{\link{gini}} for details on the components.
#'
#' @note This function gives reasonable variance estimates for basic sample
#' designs such as simple random sampling or stratified simple random sampling.
#'
#' @author Andreas Alfons
#'
#' @seealso \code{\link{variance}}, \code{\link{calibWeights}},
#' \code{\link{arpr}}, \code{\link{qsr}}, \code{\link{rmpg}}, \code{\link{gini}}
#'
#' @references
#' A. Alfons and M. Templ (2013) Estimation of Social Exclusion Indicators
#' from Complex Surveys: The \R Package \pkg{laeken}.  \emph{Journal of
#' Statistical Software}, \bold{54}(15), 1--25.  URL
#' \url{http://www.jstatsoft.org/v54/i15/}
#'
#' @keywords survey
#'
#' @examples
#' data(eusilc)
#' a <- arpr("eqIncome", weights = "rb050", data = eusilc)
#'
#' ## naive bootstrap
#' bootVar("eqIncome", weights = "rb050", design = "db040",
#'     data = eusilc, indicator = a, R = 50,
#'     bootType = "naive", seed = 123)
#'
#' ## bootstrap with calibration
#' bootVar("eqIncome", weights = "rb050", design = "db040",
#'     data = eusilc, indicator = a, R = 50,
#'     X = calibVars(eusilc$db040), seed = 123)
#'
#' @export
#' @import boot

bootVar <- function(inc, weights = NULL, years = NULL, breakdown = NULL,
                    design = NULL, cluster = NULL, data = NULL, indicator,
                    R = 100, bootType = c("calibrate", "naive"), X,
                    totals = NULL, ciType = c("perc", "norm", "basic"),
                    # type "stud" and "bca" are currently not allowed
                    alpha = 0.05, seed = NULL, na.rm = FALSE, gender = NULL,
                    method = NULL, ...) {
  UseMethod("bootVar", indicator)
}


## class "indicator"
#' @export
bootVar.indicator <- function(inc, weights = NULL, years = NULL,
                              breakdown = NULL, design = NULL, cluster = NULL,
                              data = NULL, indicator, R = 100,
                              bootType = c("calibrate", "naive"), X,
                              totals = NULL,
                              ciType = c("perc", "norm", "basic"),
                              # type "stud" and "bca" are currently not allowed
                              alpha = 0.05, seed = NULL, na.rm = FALSE,
                              gender = NULL, method = NULL, ...) {
  ## initializations
  # check whether weights have been supplied
  haveWeights <- !is.null(weights)
  haveGender <- !is.null(gender)
  # check whether indicator is broken down by year
  # if so, check whether years have been supplied
  ys <- indicator$years
  byYear <- !is.null(ys)
  if(byYear && is.null(years)) stop("'years' must be supplied")
  # check whether indicator is broken down by stratum
  # if so, check whether breakdown has been supplied
  rs <- indicator$strata
  byStratum <- !is.null(rs)
  if(byStratum && is.null(breakdown)) stop("'breakdown' must be supplied")
  haveDesign <- !is.null(design)
  haveCluster <- !is.null(cluster)
  # if a data.frame has been supplied, extract the respective vectors
  if(!is.null(data)) {
    inc <- data[, inc]
    # make numeric if indicator is proportion
    inc <- as.numeric(as.integer(inc))
    if(!is.null(weights)) weights <- data[, weights]
    if(!is.null(gender)) gender <- data[, gender]
    if(byYear) years <- data[, years]
    if(byStratum) breakdown <- data[, breakdown]
    if(haveDesign) design <- data[, design]
    if(haveCluster) cluster <- data[, cluster]
  }
  # check whether the vectors have the correct type
  # make numeric if indicator is proportion
  inc <- as.numeric(as.integer(inc))
  if(!is.numeric(inc)) stop("'inc' must be a numeric vector")
  n <- length(inc)
  if(haveWeights && !is.numeric(weights)) {
    stop("'weights' must be a numeric vector")
  }
  #	if(haveGender && !is.numeric(gender)) {
  #		stop("'gender' must be a numeric vector")
  #	}
  if(byYear && !is.numeric(years)) {
    stop("'years' must be a numeric vector")
  }
  if(byStratum && !is.vector(breakdown) && !is.factor(breakdown)) {
    stop("'breakdown' must be a vector or factor")
  }
  if(haveDesign && !is.integer(design) && !is.factor(design)) {
    stop("'design' must be an integer vector or factor")
  }
  if(haveCluster && !is.integer(cluster) && !is.factor(cluster)) {
    stop("'cluster' must be an integer vector or factor")
  }
  if(is.null(data)) {  # check vector lengths
    if(haveWeights && length(weights) != n) {
      stop("'weights' must have length ", n)
    }
    if(byYear && length(years) != n) {
      stop("'years' must have length ", n)
    }
    if(byStratum && length(breakdown) != n) {
      stop("'breakdown' must have length ", n)
    }
    if(haveDesign && length(design) != n) {
      stop("'design' must have length ", n)
    }
    if(haveCluster && length(cluster) != n) {
      stop("'cluster' must have length ", n)
    }
  }
  if(!haveDesign) design <- rep.int(1, n)
  # check other input
  if(!is.numeric(R) || length(R) == 0) stop("'R' must be numeric")
  else R <- as.integer(R[1])
  if(!is.numeric(alpha) || length(alpha) == 0) stop("'alpha' must be numeric")
  else alpha <- alpha[1]
  bootType <- match.arg(bootType)
  calibrate <- haveWeights && bootType == "calibrate"
  if(calibrate) {
    X <- as.matrix(X)
    if(!is.numeric(X)) stop("'X' must be a numeric matrix")
    # 		if(nrow(X) != n) stop("'X' must have ", n, " rows")
    if(is.null(totals)) {
      # compute totals from original data with Horvitz-Thompson estimator
      if(byYear) {
        totals <- lapply(ys,
                         function(y) {
                           # extract current year from calibration variables and
                           # weights
                           i <- years == y
                           X <- X[i, , drop=FALSE]
                           weights <- weights[i]
                           # compute totals for current year
                           apply(X, 2, function(i) sum(i*weights))
                         })
        totals <- do.call(rbind, totals)  # form matrix of totals
        rownames(totals) <- ys  # use years as rownames for totals
      } else totals <- apply(X, 2, function(i) sum(i*weights))
    } else if(byYear) totals <- as.matrix(totals)
    if(!is.numeric(totals)) stop("'totals' must be of type numeric")
  } else {
    X <- NULL
    totals <- NULL
  }
  ciType <- match.arg(ciType)

  ## preparations
  data <- data.frame(inc=inc)
  data$weight <- weights
  data$year <- years
  data$stratum <- breakdown
  data$cluster <- cluster
  data$gender <- gender
  data$method <- method  # this is a bit of an ugly hack
  if(inherits(indicator, "arpr")) {
    p <- indicator$p  # percentage of median used for threshold
  } else p <- NULL
  byP <- length(p) > 1
  if(!is.null(seed)) set.seed(seed)  # set seed of random number generator
  if(!exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)) runif(1)
  seed <- .Random.seed  # seed is later on added to the object to be returned

  ## calculations
  # get basic function for bootstrap replications with definition:
  # function(x, i, p, X, totals, rs, na.rm)
  fun <- getFun(indicator, byStratum)
  bootFun <- getBootFun(calibrate, fun)
  if(byYear) {
    # ---------- breakdown by year ----------
    # get more complex function for additional with definition
    # function(y, x, R, p, aux, totals, rs, alpha, ciType, na.rm, ...)
    funByYear <- getFunByYear(byStratum, calibrate, bootFun)
    if(byStratum) {
      # ---------- breakdown by stratum ----------
      tmp <- lapply(ys, funByYear, data, R, design, cluster, p, X, totals,
                    ys, rs, alpha, ciType, na.rm, ...)
      var <- do.call(c, lapply(tmp, function(x) x[[1]]))
      names(var) <- ys
      varByStratum <- do.call(rbind, lapply(tmp, function(x) x[[2]]))
      ci <- do.call(rbind, lapply(tmp, function(x) x[[3]]))
      rownames(ci) <- ys
      ciByStratum <- do.call(rbind, lapply(tmp, function(x) x[[4]]))
      # order 'varByStratum' and 'ciByStratum' according to 'valueByStratum'
      tmp <- indicator$valueByStratum[, 1:2]
      tmp <- data.frame(order=1:nrow(tmp), tmp)
      varByStratum <- merge(varByStratum, tmp, all=TRUE, sort=FALSE)
      varByStratum <- varByStratum[order(varByStratum$order), -4]
      ciByStratum <- merge(ciByStratum, tmp, all=TRUE, sort=FALSE)
      ciByStratum <- ciByStratum[order(ciByStratum$order), -5]
    } else {
      # ---------- no breakdown by stratum ----------
      tmp <- sapply(ys, funByYear, data, R, design, cluster, p, X, totals,
                    ys, rs, alpha, ciType, na.rm, ...)
      colnames(tmp) <- ys
      var <- tmp[1,]
      ci <- t(tmp[2:3,])
    }
  } else {
    # ---------- no breakdown by year or threshold ----------
    b <- clusterBoot(data, bootFun, R, strata=design, cluster=cluster, p=p,
                     aux=X, totals=totals, rs=rs, na.rm=na.rm, ...)
    if(byStratum) {
      # ---------- breakdown by stratum ----------
      var <- apply(b$t, 2, var)
      ci <- lapply(1:length(b$t0),
                   function(i) {
                     ci <- boot.ci(b, conf=1-alpha, type=ciType, index=i)
                     switch(ciType, perc=ci$percent[4:5], norm=ci$normal[2:3],
                            basic=ci$basic[4:5], stud=ci$student[4:5],
                            bca=ci$bca[4:5])
                   })
      ci <- do.call(rbind, ci)
      colnames(ci) <- c("lower", "upper")
      if(byP) {
        overall <- 1:length(p)
        tmp <- indicator$valueByStratum[, 1:2]
        varByStratum <- data.frame(tmp, var=var[-overall])
        var <- var[overall]
        ciByStratum <- data.frame(tmp, ci[-overall, , drop=FALSE])
        ci <- ci[overall, , drop=FALSE]
        names(var) <- rownames(ci) <- names(indicator$value)
      } else {
        varByStratum <- data.frame(stratum=rs, var=var[-1])
        var <- var[1]
        ciByStratum <- data.frame(stratum=rs, ci[-1, , drop=FALSE])
        ci <- ci[1,]
      }
    } else {
      # ---------- no breakdown by stratum ----------
      if(byP) {
        var <- apply(b$t, 2, var)
        ci <- lapply(1:length(b$t0),
                     function(i) {
                       ci <- boot.ci(b, conf=1-alpha, type=ciType, index=i)
                       switch(ciType, perc=ci$percent[4:5], norm=ci$normal[2:3],
                              basic=ci$basic[4:5], stud=ci$student[4:5],
                              bca=ci$bca[4:5])
                     })
        ci <- do.call(rbind, ci)
        colnames(ci) <- c("lower", "upper")
        names(var) <- rownames(ci) <- names(indicator$value)
      } else {
        var <- var(b$t[, 1])
        ci <- boot.ci(b, conf=1-alpha, type=ciType)
        ci <- switch(ciType, perc=ci$percent[4:5], norm=ci$normal[2:3],
                     basic=ci$basic[4:5], stud=ci$student[4:5], bca=ci$bca[4:5])
        names(ci) <- c("lower", "upper")
      }
    }
  }

  ## modify and return object
  indicator$varMethod <- "bootstrap"
  indicator$var <- var
  indicator$ci <- ci
  if(byStratum) {
    indicator$varByStratum <- varByStratum
    indicator$ciByStratum <- ciByStratum
  }
  indicator$alpha <- alpha
  indicator$seed <- seed
  return(indicator)
}


## function to perform clustered bootstrap sampling

clusterBoot <- function(data, statistic, ..., strata, cluster = NULL) {
  if(is.null(cluster)) boot(data, statistic, ..., strata=strata)
  else {
    fun <- function(cluster, i, ..., .data, .statistic) {
      # retrieve sampled individuals
      i <- do.call(c, split(1:nrow(.data), .data$cluster)[i])
      # call the original statistic for the sample of individuals
      .statistic(.data, i, ...)
    }
    keep <- !duplicated(cluster)
    boot(cluster[keep], fun, ..., strata=strata[keep],
         .data=data, .statistic=statistic)
  }
}


## utility functions: return functions to be used in the bootstrap replications


# basic function for breakdown by stratum
getFun <- function(indicator, byStratum) UseMethod("getFun")

getFun.arpr <- function(indicator, byStratum) {
  if(byStratum) {
    function(x, p, rs, na.rm) {
      threshold <- p * weightedMedian(x$inc, x$weight)
      value <- weightedRate(x$inc, x$weight, threshold, na.rm=na.rm)
      valueByStratum <- sapply(rs, function(r, x, t) {
        i <- x$stratum == r
        weightedRate(x$inc[i], x$weight[i], t, na.rm=na.rm)
      }, x=x, t=threshold)
      c(value, valueByStratum)
    }
  } else {
    function(x, p, rs, na.rm) {
      threshold <- p * weightedMedian(x$inc, x$weight)
      weightedRate(x$inc, x$weight, threshold, na.rm=na.rm)
    }
  }
}

# the argument 'p' is not necessary here, but is used so
# that we have a unified function call for all indicators
getFun.qsr <- function(indicator, byStratum) {
  if(byStratum) {
    function(x, p, rs, na.rm) {
      value <- quintileRatio(x$inc, x$weight, na.rm=na.rm)
      valueByStratum <- sapply(rs, function(r, x, t) {
        i <- x$stratum == r
        quintileRatio(x$inc[i], x$weight[i], na.rm=na.rm)
      }, x=x)
      c(value, valueByStratum)
    }
  } else {
    function(x, p, rs, na.rm) {
      quintileRatio(x$inc, x$weight, na.rm=na.rm)
    }
  }
}

getFun.rmpg <- function(indicator, byStratum) {
  if(byStratum) {
    function(x, p, rs, na.rm) {
      threshold <- 0.6 * weightedMedian(x$inc, x$weight)
      value <- relativeGap(x$inc, x$weight,
                           threshold=threshold, na.rm=na.rm)
      valueByStratum <- sapply(rs, function(r, x, t) {
        i <- x$stratum == r
        relativeGap(x$inc[i], x$weight[i], threshold=t, na.rm=na.rm)
      }, x=x, t=threshold)
      c(value, valueByStratum)
    }
  } else {
    function(x, p, rs, na.rm) {
      threshold <- 0.6 * weightedMedian(x$inc, x$weight)
      relativeGap(x$inc, x$weight, threshold=threshold, na.rm=na.rm)
    }
  }
}

# the argument 'p' is not necessary here, but is used so
# that we have a unified function call for all indicators
getFun.gini <- function(indicator, byStratum) {
  if(byStratum) {
    function(x, p, rs, na.rm) {
      value <- giniCoeff(x$inc, x$weight, na.rm=na.rm)
      valueByStratum <- sapply(rs, function(r, x, t) {
        i <- x$stratum == r
        giniCoeff(x$inc[i], x$weight[i], na.rm=na.rm)
      }, x=x)
      c(value, valueByStratum)
    }
  } else {
    function(x, p, rs, na.rm) {
      giniCoeff(x$inc, x$weight, na.rm=na.rm)
    }
  }
}

# the argument 'p' is not necessary here, but is used so
# that we have a unified function call for all indicators
getFun.prop <- function(indicator, byStratum) {
  if(byStratum) {
    function(x, p, rs, na.rm) {
      value <- propCoeff(x$inc, x$weight, na.rm=na.rm)
      valueByStratum <- sapply(rs, function(r, x, t) {
        i <- x$stratum == r
        propCoeff(x$inc[i], x$weight[i], na.rm=na.rm)
      }, x=x)
      c(value, valueByStratum)
    }
  } else {
    function(x, p, rs, na.rm) {
      propCoeff(x$inc, x$weight, na.rm=na.rm)
    }
  }
}

# the argument 'p' is not necessary here, but is used so
# that we have a unified function call for all indicators
getFun.gpg <- function(indicator, byStratum) {
  if(byStratum) {
    function(x, p, rs, na.rm) {
      value <- genderGap(x$inc, x$gender, x$method[1], x$weight, na.rm=na.rm)
      valueByStratum <- sapply(rs, function(r, x, t) {
        i <- x$stratum == r
        genderGap(x$inc[i], x$gender[i], x$method[1], x$weight[i], na.rm=na.rm)
      }, x=x)
      c(value, valueByStratum)
    }
  } else {
    function(x, p, rs, na.rm) {
      genderGap(x$inc, x$gender, x$method[1], x$weight, na.rm=na.rm)
    }
  }
}


# function that incorporates resampling and (if requested) calibration
getBootFun <- function(calibrate, fun) {
  if(calibrate) {
    function(x, i, p, aux, totals, rs, na.rm, ...) {
      x <- x[i, , drop=FALSE]
      aux <- aux[i, , drop=FALSE]
      g <- calibWeights(aux, x$weight, totals, ...)
      x$weight <- g * x$weight
      fun(x, p, rs, na.rm)
    }
  } else {
    function(x, i, p, aux, totals, rs, na.rm, ...) {
      x <- x[i, , drop=FALSE]
      fun(x, p, rs, na.rm)
    }
  }
}


# more complex function for additional breakdown by year
getFunByYear <- function(byStratum, calibrate, fun) {
  if(byStratum) {
    if(calibrate) {
      # ---------- breakdown by stratum, calibration ----------
      function(y, x, R, design, cluster, p, aux, totals, ys, rs,
               alpha, ciType, na.rm, ...) {
        i <- x$year == y
        x <- x[i, , drop=FALSE]
        aux <- aux[i, , drop=FALSE]
        design <- design[i]
        cluster <- cluster[i]
        totals <- totals[ys == y,]
        b <- clusterBoot(x, fun, R, strata=design, cluster=cluster, p=p,
                         aux=aux, totals=totals, rs=rs, na.rm=na.rm, ...)
        var <- apply(b$t, 2, var)
        varByStratum <- data.frame(year=y, stratum=rs, var=var[-1])
        var <- var[1]
        ci <- sapply(1:((length(rs) + 1)),
                     function(i) {
                       ci <- boot.ci(b, conf=1-alpha, type=ciType, index=i)
                       switch(ciType, perc=ci$percent[4:5], norm=ci$normal[2:3],
                              basic=ci$basic[4:5], stud=ci$student[4:5],
                              bca=ci$bca[4:5])
                     })
        rownames(ci) <- c("lower", "upper")
        ciByStratum <- data.frame(year=y, stratum=rs, t(ci[, -1]))
        ci <- ci[, 1]
        list(var, varByStratum, ci, ciByStratum)
      }
    } else {
      # ---------- breakdown by stratum, no calibration ----------
      function(y, x, R, design, cluster, p, aux, totals, ys, rs,
               alpha, ciType, na.rm, ...) {
        i <- x$year == y
        x <- x[i, , drop=FALSE]
        design <- design[i]
        cluster <- cluster[i]
        b <- clusterBoot(x, fun, R, strata=design, cluster=cluster, p=p,
                         aux=aux, totals=totals, rs=rs, na.rm=na.rm, ...)
        var <- apply(b$t, 2, var)
        varByStratum <- data.frame(year=y, stratum=rs, var=var[-1])
        var <- var[1]
        ci <- sapply(1:((length(rs) + 1)),
                     function(i) {
                       ci <- boot.ci(b, conf=1-alpha, type=ciType, index=i)
                       switch(ciType, perc=ci$percent[4:5], norm=ci$normal[2:3],
                              basic=ci$basic[4:5], stud=ci$student[4:5],
                              bca=ci$bca[4:5])
                     })
        rownames(ci) <- c("lower", "upper")
        ciByStratum <- data.frame(year=y, stratum=rs, t(ci[, -1]))
        ci <- ci[, 1]
        list(var, varByStratum, ci, ciByStratum)
      }
    }
  } else {
    if(calibrate) {
      # ---------- no breakdown by stratum, calibration ----------
      function(y, x, R, design, cluster, p, aux, totals, ys, rs,
               alpha, ciType, na.rm, ...) {
        i <- x$year == y
        x <- x[i, , drop=FALSE]
        aux <- aux[i, , drop=FALSE]
        design <- design[i]
        cluster <- cluster[i]
        totals <- totals[ys == y,]
        b <- clusterBoot(x, fun, R, strata=design, cluster=cluster, p=p,
                         aux=aux, totals=totals, rs=rs, na.rm=na.rm, ...)
        var <- var(b$t[, 1])
        ci <- boot.ci(b, conf=1-alpha, type=ciType)
        ci <- switch(ciType, perc=ci$percent[4:5], norm=ci$normal[2:3],
                     basic=ci$basic[4:5], stud=ci$student[4:5], bca=ci$bca[4:5])
        names(ci) <- c("lower", "upper")
        c(var, ci)
      }
    } else {
      # ---------- no breakdown by stratum, no calibration ----------
      function(y, x, R, design, cluster, p, aux, totals, ys, rs,
               alpha, ciType, na.rm, ...) {
        i <- x$year == y
        x <- x[i, , drop=FALSE]
        design <- design[i]
        cluster <- cluster[i]
        b <- clusterBoot(x, fun, R, strata=design, cluster=cluster, p=p,
                         aux=aux, totals=totals, rs=rs, na.rm=na.rm, ...)
        var <- var(b$t[, 1])
        ci <- boot.ci(b, conf=1-alpha, type=ciType)
        ci <- switch(ciType, perc=ci$percent[4:5], norm=ci$normal[2:3],
                     basic=ci$basic[4:5], stud=ci$student[4:5], bca=ci$bca[4:5])
        names(ci) <- c("lower", "upper")
        c(var, ci)
      }
    }
  }
}
