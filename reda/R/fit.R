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


#' Fit Recurrent Events Regression Based on Counts and Rate Function
#'
#' The default model is the gamma frailty model with one piece constant
#' baseline rate function, which is equivalent to negative binomial regression
#' of the same shape and rate parameter in gamma prior. 
#' Spline and piecewise constant baseline rate function can be
#' specified and applied to model fitting instead.
#' \code{rateReg} returns the fitted model
#' through a \code{\link{rateReg-class}} object.
#'
#' @details
#' Function \code{\link{Survr}} in the formula response first checks
#' the dataset and will report an error if the dataset does not
#' fall into recurrent event data framework.
#' Subject's ID is pinpointed if its observation violates any checking rule.
#' See \code{\link{Survr}} for all the checking rules.
#'
#' Function \code{rateReg} first constructs the design matrix from
#' the specified arguments: \code{formula}, \code{data}, \code{subset},
#' \code{na.action} and \code{constrasts} before model fitting.
#' The constructed design matrix will be checked again to
#' fit the recurrent event data framework
#' if any observation with missing covariates is removed.
#'
#' The model fitting process involves minimization of negative log
#' likelihood function, which calls function \code{\link[stats]{nlm}}
#' from package \pkg{stats} internally.
#' \code{help(nlm)} for more details.
#' 
#' The argument \code{start} is an optional list
#' that allows users to specify the initial guess for
#' the parameter values for the minimization of
#' negative log likelihood function.
#' The available numeric vector elements in the list include
#' \itemize{
#'     \item \code{beta}: Coefficient(s) of covariates,
#'         set to be 0.1 by default.
#'     \item \code{theta}: Parameter of frailty random effect,
#'         set to be 0.5 by default.
#'     \item \code{alpha}: Coefficient(s) of baseline rate function,
#'         set to be 0.05 by default.
#' }
#' The argument \code{control} is an optional list
#' that allows users to control the process of minimization of
#' negative log likelihood function and to specify the boundary knots,
#' intercept for baseline rate function.
#' The available elements in the list include
#' \itemize{
#'     \item \code{gradtol}: A positive scalar giving the tolerance at
#'         which the scaled gradient is considered close enough to zero
#'         to terminate the algorithm. The default value is 1e-6.
#'     \item \code{stepmax}: A positive scalar that gives the maximum
#'         allowable scaled step length. The default value is 1e5.
#'     \item \code{steptol}: A positive scalar providing the minimum
#'         allowable relative step length. The default value is 1e-6.
#'     \item \code{iterlim}: A positive integer specifying the maximum
#'         number of iterations to be performed before
#'         the program is terminated. The default value is 1e2.
#'     \item \code{boundaryKnots}: A length-two numeric vector to specify
#'         the boundary knots for baseline rate funtion. By default,
#'         the left boundary knot is zero and the right one takes the
#'         largest censoring time from data.
#'     \item \code{intercept}: A logical value specifying whether
#'         intercept is included in spline baseline rate function.
#'         For piecewise constatn baseline (\code{df}=0), the specified
#'         value would be neglected. The default value
#'         is \code{TRUE}, i.e. the intercept is included. 
#' }
#' 
#' @param formula \code{Survr} object produced by function \code{\link{Survr}}.
#' @param df An optional nonnegative integer to specify the degree of freedom
#' of baseline rate function. If argument \code{knots} or \code{degree} are
#' specified, \code{df} will be neglected whether it is specified or not.
#' @param knots An optional numeric vector that represents all the internal
#' knots of baseline rate function.
#' The default is \code{NULL}, representing no any internal knots.
#' @param degree An optional nonnegative integer to specify the degree of
#' spline bases.
#' @param data An optional data frame, list or environment containing
#' the variables in the model.  If not found in data, the variables are taken 
#' from \code{environment(formula)}, usually the environment from which 
#' function \code{\link{rateReg}} is called.
#' @param subset An optional vector specifying a subset of observations 
#' to be used in the fitting process.
#' @param na.action A function that indicates what should the procedure
#' do if the data contains \code{NA}s.  The default is set by the 
#' na.action setting of \code{\link[base]{options}}.
#' The "factory-fresh" default is \code{\link[stats]{na.omit}}.
#' Other possible values inlcude \code{\link{na.fail}},
#' \code{\link{na.exclude}}, and \code{\link{na.pass}}.
#' \code{help(na.fail)} for details.
#' @param start An optional list of starting values for the parameters
#' to be estimated in the model.  See more in section details.
#' @param control An optional list of parameters to control the
#' maximization process of negative log likelihood function
#' and adjust the baseline rate function.
#' See more in section details.
#' @param contrasts An optional list, whose entries are values 
#' (numeric matrices or character strings naming functions) to be used 
#' as replacement values for the contrasts replacement function and 
#' whose names are the names of columns of data containing factors.
#' See \code{contrasts.arg} of \code{\link[stats]{model.matrix.default}}
#' for details.
#' @param ... Other arguments for future usage.
#' @return A \code{\link{rateReg-class}} object, whose slots include
#' \itemize{
#'     \item \code{call}: Function call of \code{rateReg}.
#'     \item \code{formula}: Formula used in the model fitting.
#'     \item \code{nObs}: Number of observations.
#'     \item \code{knots}: Internal knots specified for the baseline
#'         rate function.
#'     \item \code{boundaryKnots}: Boundary knots specified for the baseline
#'         rate function.
#'     \item \code{degree}: Degree of spline bases specified in baseline
#'         rate function.
#'     \item \code{df}: Degree of freedom of the model specified.
#'     \item \code{estimates}: Estimated coefficients of covariates and
#'         baseline rate function, and estimated rate parameter of
#'         gamma frailty variable.
#'     \item \code{control}: The control list specified for model fitting.
#'     \item \code{start}: The initial guess specified for the parameters
#'         to be estimated.
#'     \item \code{na.action}: The procedure specified to deal with
#'         missing values in the covariate.
#'     \item \code{xlevels}: A list that records the levels in
#'         each factor variable.
#'     \item \code{contrasts}: Contrasts specified and used for each
#'         factor variable.
#'     \item \code{convergCode}: \code{code} returned by function
#'         \code{\link[stats]{nlm}}, which is an integer indicating why the
#'         optimization process terminated. \code{help(nlm)} for details.
#'     \item \code{logL}: Log likelihood of the fitted model.
#'     \item \code{fisher}: Observed Fisher information matrix.
#' }
#' 
#' @references 
#' Fu, H., Luo, L., & Qu Y. (2014). Hypoglycemic Events Analysis via
#' Recurrent Time-to-Event (HEART) Models. 
#' \emph{Journal of biopharmaceutical statistics}, Epub 2014 Dec 1.
#' @examples
#' library(reda)
#' 
#' ## constant rate function
#' constFit <- rateReg(Survr(ID, time, event) ~ group + x1,
#'                     data = simuDat, subset = ID %in% 1:50)
#' 
#' ## 6 pieces' piecewise constant rate function
#' piecesFit <- rateReg(Survr(ID, time, event) ~ group + x1, 
#'                      data = simuDat, subset = ID %in% 1:50,
#'                      knots = seq(28, 140, by = 28))
#'
#' ## fit rate function with cubic spline 
#' splineFit <- rateReg(Survr(ID, time, event) ~ group + x1, 
#'                      data = simuDat, subset = ID %in% 1:50,
#'                      knots = c(56, 84, 112), degree = 3)
#'
#' ## brief summary of fitted models
#' constFit
#' piecesFit
#' splineFit
#'
#' ## more specific summary
#' summary(constFit)
#' summary(piecesFit)
#' summary(splineFit)
#'
#' ## model selection based on AIC or BIC
#' AIC(constFit, piecesFit, splineFit)
#' BIC(constFit, piecesFit, splineFit)
#'
#' ## estimated covariate coefficients
#' coef(piecesFit)
#' coef(splineFit)
#'
#' ## confidence intervals for covariate coefficients
#' confint(piecesFit)
#' confint(splineFit, "x1", 0.9)
#' confint(splineFit, 1, 0.975)
#'
#' ## estimated coefficients for baseline rate function
#' baseRate(piecesFit)
#' baseRate(splineFit)
#'
#' ## estimated baseline mean cumulative function (MCF) from a fitted model
#' piecesMcf <- mcf(piecesFit)
#' plotMcf(piecesMcf, conf.int = TRUE, col = "blueviolet") +
#'     ggplot2::xlab("Days") + ggplot2::theme_bw()
#'
#' ## estimated MCF for given new data
#' newDat <- data.frame(x1 = rep(0, 2), group = c("Treat", "Contr"))
#' splineMcf <- mcf(splineFit, newdata = newDat, groupName = "Group",
#'                  groupLevels = c("Treatment", "Control"))
#' plotMcf(splineMcf, conf.int = TRUE, lty = c(1, 5)) +
#'     ggplot2::xlab("Days") + ggplot2::theme_bw()
#' 
#' @seealso
#' \code{\link{summary,rateReg-method}} for summary of fitted model;
#' \code{\link{coef,rateReg-method}} for estimated covariate coefficients;
#' \code{\link{confint,rateReg-method}} for confidence interval of
#' covariate coefficients;
#' \code{\link{baseRate,rateReg-method}} for estimated coefficients of baseline
#' rate function;
#' \code{\link{mcf,rateReg-method}} for estimated MCF from a fitted model;
#' \code{\link{plotMcf,rateRegMcf-method}} for plotting estimated MCF.
#' @importFrom splines bs
#' @importFrom stats na.fail na.omit na.exclude na.pass .getXlevels
#' @export
rateReg <- function (formula, df = NULL, knots = NULL, degree = 0L,
                     data, subset, na.action, start = list(), control = list(),
                     contrasts = NULL, ...) {
    ## record the function call to return
    Call <- match.call()

    ## arguments check
    if (missing(formula)) {
        stop("Argument 'formula' is required.")
    } 
    if (missing(data)) {
        data <- environment(formula)
    }
    if (! with(data, inherits(eval(Call[[2]][[2]]), "Survr"))) {
        stop("Response in formula must be a survival recurrent object.")
    }

    ## Prepare data: ID, time, event ~ X(s)
    mcall <- match.call(expand.dots = FALSE)
    mmcall <- match(c("formula", "data", "subset", "na.action"),
                    names(mcall), 0L)
    mcall <- mcall[c(1L, mmcall)]
    ## drop unused levels in factors 
    mcall$drop.unused.levels <- TRUE
    mcall[[1L]] <- quote(stats::model.frame)
    mf <- eval(mcall, parent.frame())
    mt <- attr(mf, "terms")
    mm <- stats::model.matrix(formula, data = mf, contrasts.arg = contrasts)
    ## get data.frame if na.action = na.pass for further data checking 
    mcall$na.action <- na.pass
    mf_na <- eval(mcall, parent.frame())
    mm_na <- stats::model.matrix(formula, data = mf_na,
                                 contrasts.arg = contrasts)
    ## number of covariates excluding intercept
    if ((nBeta <- ncol(mm) - 1L) <= 0) {
        stop("Covariates must be specified in formula.")
    }
    ## covariates' names
    covar_names <- colnames(mm)[-1]
    ## data 
    dat <- as.data.frame(cbind(mf[, 1][, 1:3], mm[, -1]))
    colnames(dat) <- c("ID", "time", "event", covar_names)
    nObs <- nrow(dat)

    ## check the impact caused by missing value
    ## if there is missing value removed
    if (nrow(mm_na) > nObs) {
        ## recover original ID names for possible pin-point
        idFactor <- with(data, attr(eval(Call[[2]][[2]]), "ID"))
        attr(dat, "ID") <- factor(levels(idFactor)[dat$ID],
                                  levels = levels(idFactor)) 
        message("Observations with missing values on covariates are removed.") 
        message("Checking new data set again ... ", appendLF = FALSE)
        check_Survr(dat)
        message("done.")
    }

    ## 'control' for 'nlm' and 'bs'
    control <- c(control, list(time = dat$time))
    control <- do.call("rateReg_control", control)
    boundaryKnots <- control$boundaryKnots
    indIntercept <- control$intercept

    ## check and reformat 'degree' at the same time
    if ((degree <- as.integer(degree)) < 0) {
        stop("'degree' must be a nonnegative integer.")
    }

    ## generate knots if knots is unspecified
    if (degree == 0L) { ## if piece-wise constant
        templist <- pieceConst(x = dat$time,
                               df = df, knots = knots)
        knots <- templist$knots
        df <- templist$df
    } else { ## else degree > 0, call 'bs' for spline 
        bsMat <- splines::bs(x = dat$time, df = df,
                             knots = knots, degree = degree,
                             intercept = indIntercept,
                             Boundary.knots = boundaryKnots)
        ## update df, knots
        knots <- as.numeric(attr(bsMat, "knots"))
        ## set bKnots as c(knots, last_boundary_knots)
        bKnots <- c(knots, boundaryKnots[2])
        df <- degree + length(knots) + as.integer(indIntercept)
        ## generate bsMat for estimated baseline rate function and mcf
        xTime <- seq(from = min(dat$time), to = max(dat$time),
                     length.out = max(1e3, length(unique(dat$time))))
        xTime <- sort(unique(c(xTime, dat$time[all.equal(dat$event, 0)])))
        bsMat_est <- splines::bs(xTime, knots = knots, degree = degree,
                                 intercept = indIntercept,
                                 Boundary.knots = boundaryKnots)
    }

    ## set bKnots as c(knots, last_boundary_knots)
    bKnots <- c(knots, boundaryKnots[2])
    alphaName <- nameBases(bKnots = bKnots, degree = degree, df = df, 
                           leftBound = boundaryKnots[1])
    
    ## start' values for 'nlm'
    startlist <- c(start, list(nBeta = nBeta, nAlpha = df))
    start <- do.call("rateReg_start", startlist)
    ini <- do.call("c", start)
    length_par <- length(ini)

    ## log likelihood
    fit <- stats::nlm(logL_rateReg, ini, data = dat, 
                      bKnots = bKnots, boundaryKnots = boundaryKnots,
                      degree = degree, bsMat = bsMat,
                      bsMat_est = bsMat_est, xTime = xTime,
                      hessian = TRUE,
                      gradtol = control$gradtol, stepmax = control$stepmax,
                      steptol = control$steptol, iterlim = control$iterlim)

    ## estimates for beta
    est_beta <- matrix(NA, nrow = nBeta, ncol = 5)
    colnames(est_beta) <- c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)")
    rownames(est_beta) <- covar_names
    
    se_vec <- sqrt(diag(solve(fit$hessian)))
    est_beta[, 1] <- fit$estimate[1:nBeta]
    est_beta[, 2] <- exp(est_beta[, 1])
    est_beta[, 3] <- se_vec[1:nBeta]
    est_beta[, 4] <- est_beta[, 1]/est_beta[, 3]
    est_beta[, 5] <- 2 * stats::pnorm(- abs(est_beta[, 4]))

    ## estimates for theta
    est_theta <- matrix(NA, nrow = 1, ncol = 2)
    colnames(est_theta) <- c("parameter", "se")
    rownames(est_theta) <- "Frailty"
    est_theta[1, ] <- c(fit$estimate[nBeta + 1], se_vec[nBeta + 1])

    ## estimates for alpha
    est_alpha <- matrix(NA, nrow = df, ncol = 2)
    colnames(est_alpha) <- c("coef", "se(coef)")
    rownames(est_alpha) <- alphaName 
    est_alpha[, 1] <- fit$estimate[(nBeta + 2):length_par]
    est_alpha[, 2] <- se_vec[(nBeta + 2):length_par]

    ## output: na.action
    if (is.null(attr(mf, "na.action"))) {
        na.action <- options("na.action")[[1]]
    } else {
        na.action <- paste("na", class(attr(mf, "na.action")), sep = ".")
    }
    
    ## output: contrasts
    if (is.null(contrasts)) {
        contrasts <- list(contrasts = NULL)
    } else {
        contrasts <- attr(mm, "contrasts")    
    }

    ## output: df, degree of freefom, including beta and theta
    df <- list(beta = nBeta, theta = 1L, alpha = df)
    
    ## results to return
    results <- methods::new("rateReg", 
                            call = Call,
                            formula = formula,
                            nObs = nObs,
                            knots = knots,
                            boundaryKnots = boundaryKnots,
                            degree = degree,
                            df = df,
                            estimates = list(beta = est_beta, 
                                             theta = est_theta, 
                                             alpha = est_alpha),
                            control = control,
                            start = start, 
                            na.action = na.action,
                            xlevels = .getXlevels(mt, mf),
                            contrasts = contrasts,
                            convergCode = fit$code,
                            logL = - fit$minimum,
                            fisher = fit$hessian)
    ## return
    results
}


## internal functions ==========================================================
whereT <- function (tt, bKnots) {
    ## designed to be used inside function 'sapply', 'tt' is length one vector
    ## return the baseline segment number 'tt' belongs to
    ## or the row number of bsMat for censoring time
    min(which(tt <= bKnots))
}

pieceConst <- function (x, df = NULL, knots = NULL) {
    ind <- (is.null(df) + 1) * is.null(knots) + 1
    ## ind == 1: knots is not NULL; df <- length(knots) + 1
    ## ind == 2: df is not NULL, while knots is NULL; number of piece <- df
    ## ind == 3: both df and knots are NULL; one-piece constant, df <- 1
    df <- switch(ind, length(knots) + 1L, as.integer(df), 1L)
    if (ind > 1) {
        tknots <- df + 1L
        quans <- seq.int(from = 0, to = 1,
                         length.out = tknots)[-c(1L, tknots)]
        knots <- as.numeric(stats::quantile(x, quans))
    }
    ## return
    list(df = df, knots = knots)
}

## baseline rate function
rho_0 <- function (par_alpha, Tvec, bKnots, degree, bsMat) {
    ## if piecewise constant, degree == 0
    if (degree == 0) {
        indx <- sapply(Tvec, whereT, bKnots)
        return(par_alpha[indx])  # function ends
    } 
    ## else spline with degree >= 1
    ## return
    bsMat %*% par_alpha
}

## mean cumulative function
mu0 <- function (par_alpha, Tvec, bKnots, degree,
                 boundaryKnots, bsMat_est, xTime = NULL) {
    ## if piecewise constant, degree == 0
    if (degree == 0) {
        ## segement number of each subject
        indx <- sapply(Tvec, whereT, bKnots)
        BL_segments <- c(bKnots[1], diff(bKnots))
        ## The MCF at each time point  
        CumMean_Pieces <- stats::diffinv(BL_segments * par_alpha)[-1]  
        mu_tau <- CumMean_Pieces[indx] -
            (bKnots[indx] - Tvec) * par_alpha[indx]
        return(mu_tau)  # function ends
    }
    ## else spline with degree >= 1
    baseRate <- bsMat_est %*% par_alpha
    if (is.null(xTime)) { ## for function mcf
        stepTime <- diff(c(boundaryKnots[1], Tvec))
        mu_tau <- cumsum(baseRate) * stepTime
        return(mu_tau)
    }
    ## else for loglikehood
    stepTime <- xTime[2] - xTime[1]
    indx <- sapply(Tvec, whereT, bKnots = xTime)
    mu_tau <- sapply(indx, function (ind) {
        sum(baseRate[seq(ind)]) * stepTime
    })
    ## return
    mu_tau
}

dmu0_dalpha <- function (tt, bKnots, degree, bsMat_est, xTime) {
    ## if baseline rate function is piecewise constant
    if (degree == 0L) {
        indx <- min(which(tt <= bKnots))
        ## BL_segments 
        value <- diff(c(0, bKnots))
        n_pieces <- length(bKnots)
        ## if tt lies in the last segment
        if (indx == n_pieces) {
            value[n_pieces] <- ifelse(n_pieces == 1, tt, 
                                      tt - bKnots[n_pieces - 1])
        } else if (indx > 1) { ## if tt lies in one of the middle segments
            value[(indx + 1) : n_pieces] <- 0
            value[indx] <- tt - bKnots[indx - 1]
        } else { ## if tt lies in the first segment
            value[(indx + 1) : n_pieces] <- 0
            value[indx] <- tt
        }
        ## return and end the function
        return(value)    
    }
    ## else it is spline with degree > 1
    stepTime <- xTime[2] - xTime[1]
    indx <- min(which(tt <= xTime))
    derVec <- sapply(indx, function (ind) {
        colSums(bsMat_est[seq(ind), ]) * stepTime
    })
    ## return
    derVec
}

dl_dalpha_part1 <- function (par_alpha, indx, degree, bsMat) {
    ## if rate function is piecewise constant
    if (degree == 0L) {
        return(1 / par_alpha * table(indx))     
    }
    ## else rate function is spline
    drho0_dbeta_ij <- bsMat
    rho0_ij <- bsMat %*% par_alpha
    ## return
    colSums(t(1 / rho0_ij) %*% drho0_dbeta_ij)
}

## compute negative log likelihood
logL_rateReg <- function (par, data, bKnots, degree,
                          boundaryKnots, bsMat, bsMat_est, xTime) {
    ## number of covariates, possibly zero
    nBeta <- ncol(data) - 3
    ## par = \THETA in the paper
    par_theta <- par[nBeta + 1]
    par_alpha <- par[(nBeta + 2) : length(par)]
    m <- length(unique(data$ID))
    xMat <- as.matrix(data[, 4:(3 + nBeta)])
    expXBeta <- exp(xMat %*% as.matrix(par[1 : nBeta]))
    ## index for event and censoring
    ind_event <- data$event == 1
    ind_cens <- data$event == 0
    ## rate function
    rho_0_ij <- rho_0(par_alpha = par_alpha,
                      Tvec = data$time[ind_event],
                      bKnots = bKnots, degree = degree, 
                      bsMat = bsMat[ind_event, ])
    rho_i <- expXBeta[ind_event] * rho_0_ij
    rho_i[rho_i < 1e-100] <- 1e-100
    sum_log_rho_i <- sum(log(rho_i))
    ## n_ij: number of event for each subject
    ## these codes to make sure that the order will not change 
    ## if the patient ID is not ordered
    n_ij <- table(data$ID)[order(unique(data$ID))] - 1  
    ## if there is a subject with 0 event, 
    ## the sequence will not be generated for this subject
    theta_j_1 <- par_theta + sequence(n_ij) - 1  
    theta_j_1[theta_j_1 < 1e-100] <- 1e-100
    sum_log_theta_j_1 <- sum(log(theta_j_1))
    ## integral that involves censoring time tau
    ## baseline mcf
    mu0i <- mu0(par_alpha = par_alpha, Tvec = data$time[ind_cens],
                bKnots = bKnots, degree = degree, boundaryKnots = boundaryKnots,
                bsMat_est = bsMat_est, xTime = xTime)
    mui <- mu0i * expXBeta[ind_cens]
    mui_theta <- par_theta + mui
    mui_theta[mui_theta < 1e-100] <- 1e-100
    sum_log_theta_mui <- sum((n_ij + par_theta) * log(mui_theta))
    if (par_theta < 1e-100) {
        par_theta <- 1e-100
    }
    logLH <- m * par_theta * log(par_theta) + sum_log_rho_i + 
        sum_log_theta_j_1 - sum_log_theta_mui
    penal <- ifelse(par_theta < 0 | min(par_alpha) < 0, 1e+50, 0)
    negLH <- -logLH + penal
    ## Calculate the gradient
    xMat_i <- xMat[ind_cens, ]
    dl_dbeta_i <- sweep(x = as.matrix(xMat_i), MARGIN = 1, FUN = "*", 
                        STATS = (n_ij - mui)/(par_theta + mui) * par_theta)
    dl_dbeta <- colSums(dl_dbeta_i)
    dl_dtheta <- m + m * log(par_theta) + 
        sum(1/(par_theta + sequence(n_ij) - 1)) - 
        sum((n_ij + par_theta)/(par_theta + mui)) - sum(log(mui_theta))
    indx <- sapply(data$time[ind_event], whereT, bKnots)
    if (length(unique(indx)) < length(bKnots)) {
        stop("Some segements have zero events!")
    }
    ## reform dimension by 'array' for one-piece baseline 
    dim_n1 <- length(par_alpha)
    dim_n2 <- length(data$time[ind_cens])
    tempPart2 <- array(sapply(data$time[ind_cens], dmu0_dalpha,
                              bKnots, degree, bsMat_est, xTime),
                       c(dim_n1, dim_n2))
    dl_dalpha_part2 <- sweep(t(tempPart2), MARGIN = 1, FUN = "*", 
                             STATS = (n_ij + par_theta) / (par_theta + mui) *
                                 expXBeta[ind_cens])
    tempPart2 <- colSums(dl_dalpha_part2)
    tempPart1 <- dl_dalpha_part1(par_alpha, indx, degree, bsMat[ind_event, ])
    dl_dalpha <-  tempPart1 - tempPart2
    attr(negLH, "gradient") <- -c(dl_dbeta, dl_dtheta, dl_dalpha)
    ## return
    negLH
}

rateReg_control <- function (gradtol = 1e-6, stepmax = 1e5, 
                             steptol = 1e-6, iterlim = 1e2,
                             boundaryKnots = NULL, intercept = TRUE,
                             time) {
    ## controls for function stats::nlm
    if (!is.numeric(gradtol) || gradtol <= 0) {
        stop("value of 'gradtol' must be > 0")
    }
    if (!is.numeric(stepmax) || stepmax <= 0) {
        stop("value of 'stepmax' must be > 0")
    } 
    if (!is.numeric(steptol) || steptol <= 0) {
        stop("value of 'steptol' must be > 0")
    } 
    if (!is.numeric(iterlim) || iterlim <= 0) {
        stop("maximum number of iterations must be > 0")
    }
    if (is.null(boundaryKnots)) {
        boundaryKnots <- c(0, max(time))
    } else {
        boundaryKnots <- sort(boundaryKnots)
        ind1 <- boundaryKnots[1] > min(time)
        ind2 <- boundaryKnots[2] < max(time)
        if (ind1 || ind2) {
            stop("boundary knots should not lie inside the range of visit time")
        }
    }
    ## return
    list(gradtol = gradtol, stepmax = stepmax, 
         steptol = steptol, iterlim = iterlim,
         boundaryKnots = boundaryKnots, intercept = intercept)
}

rateReg_start <- function (beta, theta = 0.5, alpha, nBeta, nAlpha) {
    ## beta = starting value(s) for coefficients of covariates
    ## theta = starting value for random effects
    ## alpha = starting values for piece-wise baseline rate functions
    if (missing(beta)) {
        beta <- rep(0.1, nBeta)
    } else if (length(beta) != nBeta) {
        stop(paste("number of starting values for coefficients of covariates",
                   "does not match with the specified formula"))
    }
    if (theta <= 0) {
        stop("value of parameter for random effects must be > 0")
    }
    if (missing(alpha)) {
        alpha <- rep(0.05, nAlpha)
    }
    ## return
    list(beta = beta, theta = theta, alpha = alpha)
}

## generate intervals from specified baseline pieces
nameBases <- function (bKnots, degree, df, leftBound) {
    nAlpha <- length(bKnots)
    intervals <- rep(NA, nAlpha)
    if (degree == 0L) {
        intervals[1] <- paste0("(", leftBound, ", ", bKnots[1], "]", sep = "")
        if (nAlpha > 1) {
            for(i in 2:nAlpha){
                intervals[i] <- paste0("(", bKnots[i - 1], ", ", 
                                       bKnots[i], "]", sep = "")
            }
        }
        return(intervals)
    }
    ## else degree > 0 
    ## return
    paste("B-spline", seq(df), sep = ".")
}
