################################################################################
##
##   R package reda by Wenjie Wang, Haoda Fu, and Jun Yan
##   Copyright (C) 2015
##
##   This file is part of the R package reda.
##
##   The R package reda is free software: You can redistribute it and/or
##   modify it under the terms of the GNU General Public License as published
##   by the Free Software Foundation, either version 3 of the License, or
##   any later version (at your option). See the GNU General Public License
##   at <http://www.gnu.org/licenses/> for details.
##
##   The R package reda is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##
################################################################################


## collation after class.R
#' @include class.R 
NULL


#' Mean Cumulative Function (MCF)
#' 
#' An S4 class generic function that estimates mean cumulative function (MCF)
#' from a fitted model or computing the sample nonparametric MCF
#' (also called Nelson-Aalen estimator) from data.
#' 
#' For \code{formula} object with \code{\link{Survr}} object as response, 
#' the covariate specified at the right hand side of the formula
#' should be either 1 or any one factor variable in the data.
#' The former computes the overall sample MCF.
#' The latter computes the sample MCF for each level of
#' the factor variable specified, respectively.
#' The sample MCF is also called Nelson-Aalen nonparametric estimator
#' (Nelson, 2003) and computed on each time point from sample data.
#' The point estimate of sample MCF at each time point does not
#' assume any particular underlying model. The variance
#' of estimated MCF (ReliaWiki, 2012) at each time point is estimated
#' and the approximate confidence intervals are provided as well,
#' which are constructed based on the asymptotic normality
#' of log mean cumulative function.
#' 
#' For \code{\link{rateReg-class}} object, 
#' \code{mcf} estimates the baseline MCF and its confidence interval
#' at each time grid if argument \code{newdata} is not specified.
#' Otherwise, \code{mcf} estimates MCF and its confidence interval
#' for the given newdata based on Delta-method.
#' 
#' @param object An object used to dispatch a method.
#' @param ... Other arguments for future usage.
#' @param level An optional numeric value 
#' indicating the confidence level required. 
#' The default value is 0.95.
#' @return
#' \code{\link{sampleMcf-class}} or \code{\link{rateRegMcf-class}} object.
#' Their slots include
#' \itemize{
#'     \item \code{level}: Confidence level specified.
#'     \item \code{MCF}: Mean cumulative function at each time point.
#'     \item \code{multiGroup}: A logical value indicating whether MCF
#'         is estimated for different specified group.
#'     \item \code{newdata}: Given dataset used to estimate MCF.
#' }
#' For the meaning of other slots, see \code{\link{rateReg}}.
#' 
#' @references
#' Nelson, W. B. (2003). \emph{Recurrent events data analysis for product
#' repairs, disease recurrences, and other applications} (Vol. 10). SIAM.
#'
#' ReliaWiki. (2012, March 19). Recurrent Event Data Analysis. 
#' Retrieved November 23, 2015, from
#' \url{http://reliawiki.org/index.php/Recurrent_Event_Data_Analysis}
#' @seealso \code{\link{rateReg}} for model fitting;
#' \code{\link{plotMcf}} for plotting MCF.
#' @examples 
#' library(reda)
#'  
#' ## sample MCF
#' sampleMcf <- mcf(Survr(ID, time, event) ~ group,
#'                  data = simuDat, subset = ID %in% 1:10)
#'
#' ## plot sample MCF
#' plotMcf(sampleMcf, lty = c(1, 3), col = c("orange", "navy"),
#'         mark.time = TRUE) + ggplot2::xlab("Days") + ggplot2::theme_bw()
#'
#' ## For estimated MCF from a fitted model,
#' ## see examples given in function rateReg.
#'  
#' @export
setGeneric(name = "mcf",
           def = function(object, ...) {
               standardGeneric("mcf")
           })


#' @describeIn mcf Sample MCF from data.
#' 
#' @param data An optional data frame, list or environment containing
#' the variables in the model.  If not found in data, the variables are taken 
#' from \code{environment(formula)}, usually the environment from which 
#' the function is called.
#' @param subset An optional vector specifying a subset of observations 
#' to be used in the fitting process.
#' @param na.action A function that indicates what should the procedure do 
#' if the data contains \code{NA}s.  The default is set by the 
#' na.action setting of \code{\link[base]{options}}.
#' The "factory-fresh" default is \code{\link[stats]{na.omit}}.
#' Other possible values inlcude \code{\link{na.fail}},
#' \code{\link{na.exclude}}, and \code{\link{na.pass}}.
#' \code{help(na.fail)} for details.
#' 
#' @aliases mcf,formula-method
#' @importFrom utils tail
#' @importFrom stats na.fail na.omit na.exclude na.pass
#' @export
setMethod(f = "mcf", signature = "formula", 
          definition = function(object, data, subset, na.action, 
                                level = 0.95, ...) {
              ## record the function call 
              Call <- match.call()
              Call[[1]] <- as.name("mcf")
              names(Call) <- sub(pattern = "object", 
                                 replacement = "formula", names(Call))
              mfnames <- c("formula", "data", "subset", "na.action")
              mfind <- match(mfnames, names(Call), nomatch = 0L)
              Call$formula <- eval(object)
              Call$data <- eval(substitute(alist(data)))[[1]]
              Call$subset <- eval(substitute(alist(subset)))[[1]]
              Call$na.action <- eval(substitute(alist(na.action)))[[1]]
              mcall <- Call[c(1, mfind)]
              ## Call to return
              Call <- mcall
              
              ## drop unused levels in factors 
              mcall$drop.unused.levels <- TRUE
              ## Prepare data: ID, time, event ~ X
              mcall[[1]] <- quote(stats::model.frame)
              mm <- if(is.R()) {
                  eval.parent(mcall)
              } else {
                  eval(mcall, sys.parent())
              }
              if(missing(data)) {
                  data <- environment(object)
              }
              if (! with(data, inherits(eval(object[[2]]), "Survr"))) {
                  stop("Response in formula must be a 'Survr' object.")
              }
              Terms <- stats::terms(object)
              ord <- attr(Terms, "order")
              if (length(ord) & any(ord != 1)) {
                  stop("Interaction terms are not valid for this function.")
              }
              
              ## number of covariates excluding intercept
              nBeta <- ncol(mm) - 1L
              nSample <- nrow(mm)
              if (nBeta == 0L) { 
                  X <- covar_names <- NULL
              } else {
                  X <- mm[, 2]
                  ## covariates' names
                  covar_names <- attr(Terms, "term.labels")
              }
              ## data 
              dat <- as.data.frame(cbind(mm[, 1], X))
              colnames(dat) <- c("ID", "time", "event", covar_names)     

              ## check on level specified
              level <- level[1]
              if (level <= 0 || level >= 1) {
                  stop("Confidence level must be between 0 and 1.")
              }
              
              if (nBeta == 0L) { ## if no covariates specified
                  outDat <- sMcf(dat, level = level)
                  out <- new("sampleMcf", formula = object, 
                             MCF = outDat, multiGroup = FALSE)
                  return(out)
              }
              ## else one covariate specified
              ## argument check
              if (nBeta != 1L && ! is.factor(X)) {
                  stop("The covariate in object must be a factor or '1'.")
              } 
              ## number of levels
              num_levels <- length(levels(X))
              if (num_levels == 1) {
                  warning("The factor covariate has only one level.")
              }

              outDat <- data.frame(matrix(NA, nrow = nSample, ncol = 8))
              for (i in seq(num_levels)) {
                  subdat <- base::subset(dat, subset = X %in% levels(X)[i])
                  rowInd <- if(i == 1L) {
                      seq(nrow(subdat))
                  } else {
                      seq(from = utils::tail(rowInd, 1) + 1, by = 1,
                          length.out = nrow(subdat))                            
                  }
                  outDat[rowInd, 1:7] <- sMcf(subdat, level = level)
                  outDat[rowInd, 8] <- i
              }
              outDat[, 8] <- factor(outDat[, 8], levels = seq(num_levels),
                                    labels = levels(X))
              colnames(outDat) <- c(colnames(dat)[1:3], "MCF",
                                    "var", "lower", "upper", covar_names)

              ## output: na.action
              na.action <- if (is.null(attr(mm, "na.action"))) {
                  options("na.action")[[1]]
              } else {
                  paste("na", class(attr(mm, "na.action")), sep = ".")
              }
              
              out <- methods::new("sampleMcf",
                                  call = Call,
                                  formula = object,
                                  na.action = na.action,
                                  level = level,    
                                  MCF = outDat,
                                  multiGroup = TRUE)
              ## return
              out
          })


#' @describeIn mcf Estimated MCF from a fitted model.
#' 
#' @param newdata An optional data frame. If specified, the data frame should
#' have the same column names as the covariate names appearing in the formula
#' of original fitting.
#' @param groupName An optional length-one charactor vector to specify the
#' name for grouping each unique row in \code{newdata}, such as "gender"
#' for "male" and "female". The default value is "group".
#' @param groupLevels An optional charactor vector to specify the levels for
#' each unique row in \code{newdata}, such as "treatment" and "control".
#' The default values are capital letters starting from "A". 
#' @param control An optional list to specify the time grid 
#' where the MCF are estimated.
#' The availble elements of the control list include 
#' \code{grid}, \code{length.out}, \code{from} and \code{to}.
#' The time grid can be directly specified via element \code{grid}.
#' A dense time grid is suggested.
#' Element \code{length.out} represents the length of grid points.
#' The dafault value is 1,000.
#' Element \code{from} means the starting point of grid with default 0.
#' Element \code{to} represnts the endpoint of grid
#' with the right boundary knot as default.
#' When \code{grid} is missing, the grid will be generated 
#' by \code{seq} (from package \pkg{base})
#' with arguments \code{from}, \code{to} and \code{length.out}.
#' @aliases mcf,rateReg-method
#' @importFrom stats na.fail na.omit na.exclude na.pass
#' @export
setMethod(f = "mcf", signature = "rateReg", 
          definition = function(object, newdata, groupName, groupLevels, 
                                level = 0.95, na.action, control = list(), 
                                ...) {
              beta <- object@estimates$beta[, 1]
              alpha <- object@estimates$alpha[, 1]
              fcovnames <- as.character(object@call[[2]][[3]])
              covnames <- fcovnames[fcovnames != "+"]
              nBeta <- length(beta)
              knots <- object@knots
              degree <- object@degree
              boundaryKnots <- object@boundaryKnots
              bKnots <- c(knots, boundaryKnots[2])
              interceptInd <- object@control$intercept
              controlist <- c(control, list("bKnots" = bKnots,
                                            "boundaryKnots" = boundaryKnots))
              control <- do.call("rateReg_mcf_control", controlist)
              gridTime <- control$grid
              n_xx <- control$length.out

              if (degree > 0L) {
                  ## B-spline bases at given time grid
                  bsMat_est = bs(gridTime, degree = degree,
                                 knots = knots,
                                 Boundary.knots = boundaryKnots,
                                 intercept = interceptInd)
              }
              ## estimated baseline mcf
              estMcf <- mu0(par_alpha = alpha, Tvec = gridTime,
                            bKnots = bKnots, degree = degree,
                            boundaryKnots = boundaryKnots,
                            bsMat_est = bsMat_est)

              n_par <- nrow(object@fisher)
              ## covariance matrix of beta and alpha
              covInd <- seq(n_par)[- (nBeta + 1)]
              covPar <- solve(object@fisher)[covInd, covInd]

              ## nonsense, just to suppress Note from R CMD check --as-cran
              `(Intercept)` <- NULL
              
              ## about newdata
              tt <- stats::terms(object@formula)
              Terms <- stats::delete.response(tt)
              if (missing("na.action")) {
                  na.action <- options("na.action")[[1]]
              } 
              if (missing(newdata)) {
                  X <- matrix(rep(0, nBeta), nrow = 1)
                  colnames(X) <- covnames
                  rownames(X) <- "1"
              } else {
                  mf <- stats::model.frame(Terms, newdata,
                                           na.action = na.action, 
                                           xlev = object@xlevels)
                  if (is.null(attr(mf, "na.action"))) {
                      na.action <- options("na.action")[[1]]
                  } else {
                      na.action <- paste("na", class(attr(mf, "na.action")), 
                                         sep = ".")
                  }
                  X <- stats::model.matrix(Terms, mf,
                                           contrasts.arg =
                                               object@contrasts$constracts)
                  ## remove intercept and deplicated rows
                  X <- unique(base::subset(X, select = -`(Intercept)`))
                  if (ncol(X) != nBeta) {
                      stop(paste("The number of input covariates does not",
                                 "match with 'rateReg' object"))
                  }
              }

              ## mcf for possible multigroups
              ndesign <- nrow(X)
              multiGroup <- ifelse(ndesign == 1, FALSE, TRUE)
              
              coveff <- as.numeric(exp(X %*% beta))
              outDat <- matrix(NA, ncol = 4, nrow = ndesign * n_xx)

              for (i in seq(ndesign)) {
                  ## Delta-method
                  gradMat <- gradDelta(Xi = X[i, ], estMcf = estMcf,
                                       degree = degree, bKnots = bKnots, 
                                       boundaryKnots = boundaryKnots,
                                       gridTime = gridTime,
                                       bsMat_est = bsMat_est,
                                       coveffi = coveff[i])
                  varTime <- apply(gradMat, 1, function (gradVec, covMat) {
                      crossprod(gradVec, covMat) %*% gradVec
                  }, covMat = covPar)
                  confBand <- stats::qnorm((1 + level) / 2) * sqrt(varTime)
                  mcf_i <- estMcf * coveff[i]
                  lower <- pmax(0, mcf_i - confBand)
                  upper <- mcf_i + confBand
                  ind <- seq(from = n_xx * (i - 1) + 1,
                             to = n_xx * i, by = 1)
                  outDat[ind, 1] <- gridTime
                  outDat[ind, 2] <- mcf_i
                  outDat[ind, 3] <- lower
                  outDat[ind, 4] <- upper
              }
              outDat <- as.data.frame(outDat)
              colnames(outDat) <- c("time", "MCF", "lower", "upper")
              if (multiGroup) {
                  if (missing(groupLevels)) {
                      groupLevels <- LETTERS[seq(ndesign)]
                  }
                  if (missing(groupName)) {
                      groupName <- "group"
                  }
                  tempcol <- factor(rep(groupLevels[seq(ndesign)], each = n_xx))
                  outDat <- cbind(outDat, tempcol)
                  colnames(outDat) <- c("time", "MCF", "lower", "upper",
                                        groupName)
              }
              ## output
              out <- new("rateRegMcf",
                         call = object@call,
                         formula = object@formula,
                         knots = knots, degree = degree, 
                         boundaryKnots = boundaryKnots,
                         newdata = X, MCF = outDat, level = level,
                         na.action = na.action, control = control, 
                         multiGroup = multiGroup)
              ## return
              out
          })


## internal function ===========================================================
rateReg_mcf_control <- function (grid, length.out = 1e3, from, to, 
                                 bKnots, boundaryKnots) {
    ## controls for function MCF with signiture rateReg
    from <- if (missing(from)) boundaryKnots[1]
    to  <- if (missing(to)) boundaryKnots[2]
    if (! missing(grid)) {
        if (! is.numeric(grid) || is.unsorted(grid)) {
            stop("'grid' specified must be a increasing numeric vector.")
        }
        length.out <- length(grid)
        from <- min(grid)
        to <- max(grid)  
    } else {
        grid <- seq(from = from, to = to, length.out = length.out)
    }
    if (min(grid) < boundaryKnots[1] || max(grid) > boundaryKnots[2]) {
        stop("'grid' must be within the coverage of boundary knots.")
    }
    ## return
    list(grid = grid, length.out = length.out, from = from, to = to)
}

## compute sample MCF
sMcf <- function(inpDat, level) {
    ## if time ties, put event time before censoring time
    tempDat <- with(inpDat, inpDat[base::order(time, 1 - event), 1:3])
    nSample <- nrow(tempDat)
    num_pat <- length(unique(tempDat$ID))
    num_at_risk <- num_pat - cumsum(tempDat$event == 0)
    increment <- ifelse(tempDat$event == 1, 1 / num_at_risk, 0)
    ## index of unique event and censoring time
    indx <- ! duplicated(tempDat$time, fromLast = TRUE)
    ## sample mcf at each time point
    smcf <- cumsum(increment)
    inc2 <- increment ^ 2
    inc12 <- (1 - increment) ^ 2
    incre <- inc2 * (inc12 + (num_at_risk - 1) * inc2)
    incre_var <- ifelse(tempDat$event == 0, 0, incre)
    var_smcf <- cumsum(incre_var)
    wtonexp <- stats::qnorm(0.5 + level / 2) * sqrt(var_smcf) / smcf
    upper <- smcf * exp(wtonexp)
    lower <- smcf * exp(- wtonexp)
    ## return
    data.frame(tempDat, MCF = smcf, var = var_smcf, lower, upper)
}

## Delta-method, compute the gradient on beta and alpha
gradDelta <- function (Xi, estMcf, degree, bKnots, boundaryKnots,
                       gridTime, bsMat_est, coveffi) {
    gradBeta <- estMcf %o% Xi * coveffi
    allKnots <- c(boundaryKnots[1], bKnots)
    n_xx <- length(gridTime)
    stepTime <- diff(c(boundaryKnots[1], gridTime))
    if (degree == 0L) { ## if piecewise constant rate function
        gradAlpha <- matrix(NA, nrow = n_xx, ncol = length(bKnots))
        for (j in seq_along(bKnots)) {
            gradAlpha[, j] <- sapply(gridTime, function (tt, knotj, knotjm1) {
                (tt >= knotj) * (knotj - knotjm1) +
                    (tt < knotj && tt > knotjm1) * (tt - knotjm1)          
            }, knotj = allKnots[j + 1], knotjm1 = allKnots[j])
        }
    } else { ## degree >= 1, spline rate function
        gradAlpha <- apply(bsMat_est, 2, cumsum) * stepTime 
    }
    gradAlpha <- gradAlpha * coveffi
    ## retrun
    cbind(gradBeta, gradAlpha)
}
