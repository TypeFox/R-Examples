#' Methods for expanded datasets
#'
#' @description Regression weights, residuals and residual plots for expanded datasets.
#' 
#' @param object an expanded dataset (of class \code{"\link{expData}"}).
#' @param ... additional arguments.
#' @details 
#' \code{weights} extracts regression weights (to be used in the natural effect model) for each observation of an expanded dataset.
#' 
#' \code{residuals} extracts residuals from the working model which is stored as an attribute of the expanded dataset. These can be used to assess normality of the residuals of the mediator working model when using the weighting-based approach (see example).
#' 
#' \code{residualPlot} and \code{residualPlots} are convenience functions from the \pkg{car} package. These can be used to assess the adequacy of the working model.
#'
# The regression weights are a multiplication of (and hence reflect)
# \enumerate{
#  \item ratio-of-mediator probability (density) weights: only when the weighted-based approach is used and \code{object} hence inherits from class \code{"\link{weightData}"}
#  \item survey weights
# }
#' @name expData-methods
#' @seealso \code{\link{expData}}, \code{\link{neWeight}}, \code{\link[car]{residualPlot}}, \code{\link[car]{residualPlots}}, \code{\link{residuals}}, \code{\link{weights}}
#' @examples
#' data(UPBdata)
#' 
#' weightData <- neWeight(negaff ~ att + gender + educ + age, 
#'                        data = UPBdata, nRep = 2)
#' 
#' ## extract regression weights for natural effect model
#' head(weights(weightData)) 
#' 
#' ## assess normality
#' qqnorm(residuals(weightData))
#' 
#' ## assess model adequacy
#' library(car)
#' residualPlots(weightData)
NULL

#' Methods for natural effect models
#'
#' @description Extractor functions, confidence intervals, residual plots and statistical tests for natural effect models.
#' @param object a fitted natural effect model object.
#' @param ... additional arguments.
# (see \code{\link[boot]{boot.ci}} for \code{confint} or \code{\link[stats]{summary.glm}} for \code{summary}).
#' @inheritParams stats::confint.default
#' @inheritParams neLht-methods
#' @details 
#' \code{confint} yields bootstrap confidence intervals or confidence intervals based on the sandwich estimator (depending on the type of standard errors requested when fitting the \code{\link{neModel}} object). 
#' Bootstrap confidence intervals are internally called via the \code{\link[boot]{boot.ci}} function from the \pkg{boot} package.
#' Confidence intervals based on the sandwich estimator are internally called via \code{\link[stats]{confint.default}}.
#' The default confidence level specified in \code{level} (which corresponds to the \code{conf} argument in \code{\link[boot]{boot.ci}}) is 0.95
#' and the default type of bootstrap confidence interval, \code{"norm"}, is based on the normal approximation.
#' Bias-corrected and accelerated (\code{"bca"}) bootstrap confidence intervals require a sufficiently large number of bootstrap replicates (for more details see \code{\link[boot]{boot.ci}}).
#'
#' A summary table with large sample tests, similar to that for \code{\link[stats]{glm}} output, can be obtained using \code{summary}.
#'
#' \code{vcov} returns either the bootstrap variance-covariance matrix (calculated from the bootstrap samples stored in\cr \code{object$bootRes}; see \code{\link{neModel}})
#' or the robust variance-covariance matrix (which is a diagonal block matrix of the original sandwich covariance matrix).
#'
#' \code{weights} returns a vector containing the regression weights used to fit the natural effect model.
#' These weights can be based on
#' \enumerate{
#'  \item ratio-of-mediator probability (density) weights (only if the weighting-based approach is used)
#'  \item inverse probability of treatment (exposure) weights (only if \code{xFit} was specified in \code{\link{neModel}})
#' }
#'
#' \code{residualPlot} and \code{residualPlots} are convenience functions from the \pkg{car} package. These can be used to assess model adequacy.
#'
#' @name neModel-methods
#' @note For the bootstrap, \emph{z}-values in the summary table are calculated by dividing the parameter estimate by its corresponding bootstrap standard error. 
#' Corresponding \emph{p}-values in the summary table are indicative, since the null distribution for each statistic is assumed to be approximately standard normal.
#' Therefore, whenever possible, it is recommended to focus mainly on bootstrap confidence intervals for inference, rather than the provided \emph{p}-values.
#' @seealso \code{\link{neModel}}, \code{\link{plot.neModel}}, \code{\link[car]{residualPlot}}, \code{\link[car]{residualPlots}}, \code{\link{weights}}
#' @examples
#' data(UPBdata)
#' 
# impData <- neImpute(UPB ~ att * negaff + educ + gender + age, 
#                     family = binomial, data = UPBdata)
#' weightData <- neWeight(negaff ~ att + educ + gender + age,
#'                        data = UPBdata)
#' neMod <- neModel(UPB ~ att0 * att1 + educ + gender + age, 
#'                  family = binomial, expData = weightData, se = "robust")
#' 
#' ## extract coefficients
#' coef(neMod)
#' 
#' ## extract variance-covariance matrix
#' vcov(neMod)
#' 
#' ## extract regression weights
#' w <- weights(neMod)
#' head(w)
#' 
#' ## obtain bootstrap confidence intervals
#' confint(neMod)
#' confint(neMod, parm = c("att0"))
#' confint(neMod, type = "perc", level = 0.90)
#' 
#' ## summary table
#' summary(neMod)
#' 
#' ## residual plots
#' library(car)
#' residualPlots(neMod)
NULL

#' Confidence interval plots for natural effect components
#'
#' @description Obtain effect decomposition confidence interval plots for natural effect models.
#' @param x a fitted natural effect model object.
#' @inheritParams plot.neLht
#' @inheritParams neEffdecomp
#' @details This function yields confidence interval plots for the natural effect components.
#' These causal parameter estimates are first internally extracted from the \code{neModel} object by applying the effect decomposition function \code{\link{neEffdecomp}(x, xRef, covLev)}.
#' @name plot.neModel
#' @examples
#' data(UPBdata)
#' 
#' impData <- neImpute(UPB ~ att * negaff + educ + gender + age, 
#'                     family = binomial, data = UPBdata)
#' neMod <- neModel(UPB ~ att0 * att1 + educ + gender + age, 
#'                  family = binomial, expData = impData, se = "robust")
#'
#' plot(neMod)
#' plot(neMod, transf = exp, 
#'      ylabels = c("PDE", "TDE", "PIE", "TIE", "TE"))
#' plot(neMod, level = 0.9, xRef = c(-1, 0))
NULL

#' @rdname neModel-methods
#' @export
coef.neModel <- function (object, ...) 
{
    coef(object$neModelFit)
}

#' @rdname neModel-methods
#' @export
confint.neModelBoot <- function (object, parm, level = 0.95, type = "norm", ...) 
{
    if (missing(parm)) 
        parm <- names(coef(object))
    else if (is.numeric(parm)) 
        parm <- names(coef(object))[parm]
    extractBootci <- function(index, x, level, type) rev(rev(boot::boot.ci(x, 
        conf = level, type = type, index = index)[[4]])[1:2])
    ci <- t(sapply(1:length(object$bootRes$t0), extractBootci, 
        object$bootRes, level = level, type = type))
    dimnames(ci) <- list(names(coef(object)), paste0(100 * level, 
        c("% LCL", "% UCL")))
    ci <- ci[parm, ]
    attributes(ci) <- c(attributes(ci), list(level = level, coef = coef(object)[parm], 
        R = object$bootRes$R, type = type))
    class(ci) <- c("neBootCI", "neModelBootCI", class(ci))
    return(ci)
}

#' @rdname neModel-methods
#' @export
confint.neModel <- function (object, parm, level = 0.95, ...) 
{
    ci <- confint.default(object, parm, level, ...)
    dimnames(ci)[[2]] <- paste0(100 * level, c("% LCL", "% UCL"))
    attributes(ci) <- c(attributes(ci), list(level = level, coef = coef(object)[parm]))
    class(ci) <- c("neModelCI", class(ci))
    return(ci)
}

#' @export
df.residual.neModel <- function (object, ...) 
{
    df.residual(object$neModelFit, ...)
}

#' @export
model.matrix.neModel <- function (object, ...) 
{
    model.matrix(object$neModelFit, ...)
}

#' Fit a natural effect model
#'
#' @description \code{neModel} is used to fit a natural effect model on the expanded dataset.
#' @param formula a \code{\link[stats]{formula}} object providing a symbolic description of the natural effect model.
#' @param expData the expanded dataset (of class \code{"\link{expData}"}).
#' @param xFit fitted model object representing a model for the exposure (used for inverse treatment (exposure) probability weighting).
#' @param se character string indicating the type of standard errors to be calculated. The default type is based on the bootstrap (see details).
#' @param nBoot number of bootstrap replicates (see \code{R} argument of \code{\link[boot]{boot}}).
#' @param parallel (only for bootstrap) The type of parallel operation to be used (if any). If missing, the default is taken from the option \code{"boot.parallel"} (and if that is not set, \code{"no"}).
#' @param ncpus (only for bootstrap) integer: number of processes to be used in parallel operation: typically one would chose this to the number of available CPUs (see details).
#' @param progress (only for bootstrap) logical value indicating whether or not a progress bar should be displayed. Progress bars are automatically disabled for multicore processing.
#' @param ... additional arguments (passed to \code{\link[stats]{glm}}).
#' @inheritParams stats::glm
# @inheritParams boot::boot
#' @return An object of class \code{"\link[=neModel-methods]{neModel}"} (which additionally inherits from class \code{"neModelBoot"} if the bootstrap is used) consisting of a list of 3 objects:
#' \item{\code{neModelFit}}{the fitted natural model object (of class \code{"\link[stats]{glm}"}) with downwardly biased standard errors}
#' \item{\code{bootRes}, \code{vcov}}{the bootstrap results (of class \code{"\link[boot]{boot}"}; if \code{se = "bootstrap"}) or the robust variance-covariance matrix (if \code{se = "robust"})}
#' \item{\code{terms}}{the \code{neTerms} (internal class) object used. This object is equivalent to the \code{\link[=terms.object]{terms}} object returned by the \code{\link[stats]{glm}} function, 
#' but has an additional \code{"vartype"} attribute, a list including pointers to the names of the outcome variable (\code{Y}), exposure (\code{X}), mediator (\code{M}), covariates (\code{C}) and auxiliary hypothetical variables \emph{x} and \emph{x*} (\code{Xexp}).}
#' See \code{\link{neModel-methods}} for methods for \code{neModel} objects.
#' 
#' @details This function is a wrapper for \code{\link[stats]{glm}}, providing unbiased bootstrap (\code{se = "bootstrap"}, the default) or robust (\code{se = "robust"}) standard errors for the parameter estimates (see below for more details).
#'
#' The \code{formula} argument requires to be specified in function of the variables from the expanded dataset (specified in \code{expData}) whose corresponding parameters index the direct and indirect effect.
#' Stratum-specific natural effects can be estimated by additionally modeling the relation between the outcome and baseline covariates.
#' If the set of baseline covariates adjusted for in the \code{formula} argument is not sufficient to control for confounding (e.g. when fitting a population-average natural effect model),
#' an adequate model for the exposure (conditioning on a sufficient set of baseline covariates) should be specified in the \code{xFit} argument.
#' In this case, such a model for the exposure distribution is needed to weight by the reciprocal of the probability (density) of the exposure (i.e. inverse probability weighting) in order to adjust for confounding.
#' Just as for ratio-of-mediator probability weighting (see paragraph below), this kind of weighting is done internally.
#'
#' Quadratic or higher-order polynomial terms can be included in the \code{formula} by making use of the \code{\link[base]{I}} function or by using the \code{\link[stats]{poly}} function.
#' However, we do not recommend the use of orthogonal polynomials (i.e. using the default argument specification \code{raw = FALSE} in \code{poly}), as these are not compatible with the \code{\link[medflex]{neEffdecomp}} function.
#'
#' In contrast to \code{\link[stats]{glm}}, the \code{expData} argument (rather than \code{data} argument) requires specification of a data frame that inherits from class \code{"\link{expData}"},
#' which contains additional information about e.g. the fitted working model, the variable types or terms of this working model 
#' and possibly ratio-of-mediator probability weights.
#' The latter are automatically extracted from the \code{\link{expData}} object and weighting is done internally.
#'
#' As the default \code{\link[stats]{glm}} standard errors fail to reflect the uncertainty inherent to the working model(s) (i.e. either a model for the mediator or an imputation model for the outcome and possibly a model for the exposure),
#' bootstrap standard errors (using the \code{\link[boot]{boot}} function from the \pkg{boot} package) or robust standard errors are calculated. The default type of standard errors is bootstrap standard errors. 
#' Robust standard errors (based on the sandwich estimator) can be requested (to be calculated) instead by specifying \code{se = "robust"}.
#' 
#' @section Bootstrap standard errors:
#' The bootstrap procedure entails refitting all working models on each bootstrap sample, reconstructing the expanded dataset and subsequently refitting the specified natural effect model on this dataset.
#' In order to obtain stable standard errors, the number of bootstrap samples (specified via the \code{nBoot} argument) should be chosen relatively high (default is 1000).
#' 
#' To speed up the bootstrap procedure, parallel processing can be used by specifying the desired type of parallel operation via the \code{parallel} argument (for more details, see \code{\link[boot]{boot}}).
#' The number of parallel processes (\code{ncpus}) is suggested to be specified explicitly (its default is 1, unless the global option \code{options("boot.cpus")} is specified). 
#' The function \code{\link[parallel]{detectCores}} from the \pkg{parallel} package can be helpful at determining the number of available cores (although this may not always correspond to the number of \emph{allowed} cores).
#'
#' @section Robust standard errors:
#' Robust variance-covariance matrices for the model parameters, based on the sandwich estimator, are calculated using core functions from the \pkg{sandwich} package.
#' Additional details and derivations for the sandwich estimator for natural effect models can be found in the corresponding vignette that can be obtained by the command \code{vignette("sandwich", package = "medflex")}.
#'
#' @note It is important to note that the original mediator(s) should not be specified in the \code{formula} argument, as the natural indirect effect in natural effect models
#' should be captured solely by parameter(s) corresponding to the auxiliary hypothetical variable \emph{x*} in the expanded dataset (see \code{\link{expData}}).
#'
#' @seealso \code{\link{neModel-methods}}, \code{\link{plot.neModel}}, \code{\link{neImpute}}, \code{\link{neWeight}}, \code{\link{neLht}}, \code{\link{neEffdecomp}}
#'
#' @examples
#' data(UPBdata)
#' 
#' ##############################
#' ## weighting-based approach ##
#' ##############################
#' weightData <- neWeight(negaff ~ att + gender + educ + age, 
#'                        data = UPBdata)
#' 
#' ## stratum-specific natural effects
#' # bootstrap SE
#' \dontrun{
#' weightFit1b <- neModel(UPB ~ att0 * att1 + gender + educ + age, 
#'                        family = binomial, expData = weightData)
#' summary(weightFit1b)
#' }
#' # robust SE
#' weightFit1r <- neModel(UPB ~ att0 * att1 + gender + educ + age, 
#'                        family = binomial, expData = weightData, se = "robust")
#' summary(weightFit1r)
#' 
#' ## population-average natural effects
#' expFit <- glm(att ~ gender + educ + age, data = UPBdata)
#' # bootstrap SE
#' \dontrun{
#' weightFit2b <- neModel(UPB ~ att0 * att1, family = binomial, 
#'                        expData = weightData, xFit = expFit)
#' summary(weightFit2b)
#' }
#' # robust SE
#' weightFit2r <- neModel(UPB ~ att0 * att1, family = binomial, 
#'                        expData = weightData, xFit = expFit, se = "robust")
#' summary(weightFit2r)
#' 
#' ###############################
#' ## imputation-based approach ##
#' ###############################
#' impData <- neImpute(UPB ~ att * negaff + gender + educ + age, 
#'                     family = binomial, data = UPBdata)
#' 
#' ## stratum-specific natural effects
#' # bootstrap SE
#' \dontrun{
#' impFit1b <- neModel(UPB ~ att0 * att1 + gender + educ + age, 
#'                     family = binomial, expData = impData)
#' summary(impFit1b)
#' }
#' # robust SE
#' impFit1r <- neModel(UPB ~ att0 * att1 + gender + educ + age, 
#'                     family = binomial, expData = impData, se = "robust")
#' summary(impFit1r)
#' 
#' ## population-average natural effects
#' # bootstrap SE
#' \dontrun{
#' impFit2b <- neModel(UPB ~ att0 * att1, family = binomial, 
#'                     expData = impData, xFit = expFit)
#' summary(impFit2b)
#' }
#' # robust SE
#' impFit2r <- neModel(UPB ~ att0 * att1, family = binomial, 
#'                     expData = impData, xFit = expFit, se = "robust")
#' summary(impFit2r)
#' 
#' \dontshow{# check with vgam (VGAM package)
#' # library(VGAM)
#' # weightData <- neWeight(negaff ~ att + gender + educ + age, family = "gaussianff", data = UPBdata, FUN = vgam)
#' # impData <- neImpute(UPB ~ att + negaff + gender + educ + age, family = "binomialff", data = UPBdata, FUN = vgam)
#' # debug(neModel)
#' # weightFit <- neModel(UPB ~ att0 + att1 + gender + educ + age, family = binomial, expData = weightData, nBoot = 2)
#' # impFit <- neModel(UPB ~ att0 + att1 + gender + educ + age, family = binomial, expData = impData, nBoot = 2)
#' # summary(weightFit)
#' # summary(impFit)
#'
#' # warning!
#' impFit <- neModel(UPB ~ att0 * att1 + gender + educ + age, family = binomial, expData = impData, nBoot = 2)
#'
#' # inverse propensity score weighting
#' expFit <- glm(att ~ gender + educ + age, data = UPBdata)
#' impFit <- neModel(UPB ~ att0 + att1, family = binomial, expData = impData, xFit = expFit, nBoot = 2)}
#' @references
#' Lange, T., Vansteelandt, S., & Bekaert, M. (2012). A Simple Unified Approach for Estimating Natural Direct and Indirect Effects. \emph{American Journal of Epidemiology}, \bold{176}(3), 190-195.
#' 
#' Vansteelandt, S., Bekaert, M., & Lange, T. (2012). Imputation Strategies for the Estimation of Natural Direct and Indirect Effects. \emph{Epidemiologic Methods}, \bold{1}(1), Article 7.
#' 
#' Loeys, T., Moerkerke, B., De Smet, O., Buysse, A., Steen, J., & Vansteelandt, S. (2013). Flexible Mediation Analysis in the Presence of Nonlinear Relations: Beyond the Mediation Formula. \emph{Multivariate Behavioral Research}, \bold{48}(6), 871-894.
#' @export
neModel <- function (formula, family = gaussian, expData, xFit, se = c("bootstrap", "robust"), 
    nBoot = 1000, parallel = c("no", "multicore", "snow"), ncpus = getOption("boot.ncpus", 
        1L), progress = TRUE, ...) 
{
    args <- as.list(match.call()) 
    if (missing(expData)) {
      expData <- parent.frame(2)$envir
      args$expData <- quote(expData)
    }
    args[[1]] <- substitute(neModelEst)
    args <- c(args[[1]], args[names(args) %in% names(formals(neModelEst))])
    neModelFit <- eval(as.call(args))
    fit <- attr(expData, "model")
    se <- match.arg(se)
    if (se == "robust") {
        progress <- FALSE
    }
    else {
        parallel <- match.arg(parallel)
        if (parallel != "no") {
          if (progress) {
            message("Progress bar disabled for parallel processing.")
            progress <- FALSE
          }
        }
    }
    if (progress) {
        env <- environment()
        counter <- 0
        progbar <- txtProgressBar(min = 0, max = nBoot, style = 3)
    }
    m <- attr(terms(expData), "vartype")$M
    mpattern <- mapply(glob2rx, c(paste0(m), paste0(m, ":*"), paste0("*:", m, ":*"), paste0("*:", m), paste0("*(", m, ")*"), paste0("*(", m, "^*")), trim.head = TRUE)
    if (sum(sapply(mpattern, grepl, glob2rx(labels(terms(neModelFit))), fixed = TRUE)))
        stop("The original mediator variables should not be included in the natural effect model!")
    if (inherits(expData, "impData")) {
        nMed <- length(attr(terms(expData), "vartype")$M)
        xasis <- attr(attr(terms(expData), "vartype"), "xasis")
        masis <- attr(attr(terms(expData), "vartype"), "masis")
        xpattern <- mapply(glob2rx, c(paste0(xasis), paste0(xasis, ":*"), paste0("*:", xasis, ":*"), paste0("*:", xasis)), trim.head = TRUE)
        mpattern <- mapply(glob2rx, c(paste0(masis), paste0(masis, ":*"), paste0("*:", masis, ":*"), paste0("*:", masis)), trim.head = TRUE)
        X0 <- attr(terms(expData), "vartype")$Xexp[1]
        X1 <- attr(terms(expData), "vartype")$Xexp[2]
        pattern <- c(xpattern, mpattern)
        termsImpData <- labels(terms(expData))
        replacement <- c(X0, paste0(X0, ":"), paste0(":", X0, ":"), paste0(":", X0), rep(c(X1, paste0(X1, ":"), paste0(":", X1, ":"), paste0(":", X1)), each = nMed))
        termsImpData <- unique(mgsub(pattern, replacement, termsImpData))
        termsNeModel <- mgsub(dimnames(attr(terms(neModelFit), 
            "factors"))[[1]], all.vars(neModelFit$formula), labels(terms(neModelFit)), 
            fixed = TRUE)
        if (sum(!termsNeModel %in% termsImpData)) 
            warning("The imputation model should at least reflect the structure of the natural effect model! This may cause some of the effect estimates to be attenuated. Please consult the terms of the models 'x' using 'labels(terms(x))'.")
    }
    switch(se, 
           "bootstrap" = {
             neModelBootfun <- function(data, ind, fit, expData, neModelFit) {
               bootData <- data[ind, ]
               FUN <- extrCall(fit)[[1]]
               call1 <- as.list(match.call(eval(FUN), extrCall(fit)))
               if (inherits(fit, "SuperLearner")) 
                 call1$X[[2]] <- call1$Y[[2]] <- substitute(bootData)
               else call1$data <- substitute(bootData)
               bootExpFit <- eval(as.call(call1))
               nRep <- nrow(expData) / nrow(data)
               indExp <- as.vector(sapply(ind, function(x) (x * nRep - 
                                                              (nRep - 1)):(x * nRep)))
               bootExpData <- expData[indExp, ]
               bootExpData <- eval(attr(expData, "call")[[1]])(bootExpFit, 
                                                               skipExpand = TRUE, expData = bootExpData)
               if (!is.null(extrCall(neModelFit)$xFit)) {
                 call3 <- as.list(extrCall(eval(extrCall(neModelFit)$xFit)))
                 call3$data <- substitute(bootData)
                 bootxFit <- eval(as.call(call3))
               }
               call4 <- as.list(attr(neModelFit, "call"))
               call4$expData <- substitute(bootExpData)
               if (!is.null(extrCall(neModelFit)$xFit)) 
                 call4$xFit <- substitute(bootxFit)
               bootFit <- eval(as.call(call4))
               if (progress) {
                 curVal <- get("counter", envir = env)
                 assign("counter", curVal + 1, envir = env)
                 setTxtProgressBar(get("progbar", envir = env), curVal + 
                                     1)
               }
               return(coef(bootFit))
             }
             bootRes <- boot::boot(data = extrData(fit), statistic = neModelBootfun, 
                                   R = nBoot, fit = fit, expData = expData, neModelFit = neModelFit, 
                                   parallel = parallel, ncpus = ncpus)
             neModel <- list(neModelFit = neModelFit, bootRes = bootRes)
           },
           "robust" = {
             fit1 <- neModelFit
             fit2 <- attr(expData, "model")
             fit3 <- if (!missing(xFit)) xFit else NULL
             
             if (inherits(expData, "weightData") && !is.null(na.action(fit1))) {
               fit1$data$wts <- vector(length = nrow(fit1$data))
               fit1$data$wts[-na.action(fit1)] <- weights(fit1, "prior")
               fit1$data[na.action(fit1), all.vars(fit1$formula)[1]] <- 0
               fit1$call[[1]] <- quote(glm)
               fit1$call[["expData"]] <- NULL
               fit1$call[["data"]] <- quote(fit1$data)
               fit1$call[["weights"]] <- quote(wts)
               if (!is.null(fit1$call[["xFit"]])) fit1$call[["xFit"]] <- NULL
               fit1 <- update(fit1)
             }
             
             if (inherits(expData, "impData") && !is.null(na.action(fit2))) {
               fit2$data$wts <- vector(length = nrow(fit2$data))
               fit2$data$wts[-na.action(fit2)] <- weights(fit2, "prior")
               fit2$data[na.action(fit2), all.vars(fit2$formula)[1]] <- 0
               fit2$call[["data"]] <- quote(fit2$data)
               fit2$call[["weights"]] <- quote(wts)
               fit2 <- update(fit2)
             }
             
             fit <- list(fit1, fit2, fit3)
             rm(fit1, fit2, fit3)
             fit <- fit[!sapply(fit, is.null)]
             
             if (any(!sapply(fit, function(x) inherits(x, "glm"))))
               stop("Robust standard errors only available if all working models are fitted via glm function. Use the bootstrap instead.")
             if (all(inherits(expData, "weightData"), !fit[[2]]$family$family %in% c("gaussian", "binomial", "poisson"))) 
               stop("Robust standard errors only available for weighting-based approach if family of the working model is either 'gaussian', 'binomial' or 'poisson'.")
             if (length(fit) > 2 && !fit[[3]]$family$family %in% c("gaussian", "binomial", "poisson"))
               stop("Robust standard errors only available if family of the exposure model is either 'gaussian', 'binomial' or 'poisson'.")

             coefnames <- lapply(fit, function(x) names(coef(x)))
             dimnames <- unlist(coefnames)
             ind <- lapply(coefnames, indmatch, dimnames)
             
             ## ESTIMATING EQUATIONS        
             estEqList <- lapply(fit, sandwich::estfun)
             estEqList[[1]] <- aggregate(estEqList[[1]], by = list(as.numeric(expData$id)), FUN = mean)[, -1]
             row.names(estEqList[[1]]) <- unique(expData$id)
             estEq <- as.matrix(data.frame(estEqList))
             rm(estEqList)
             dimnames(estEq)[[2]] <- dimnames
             
             ## MEAT
             meat <- crossprod(estEq) / nrow(estEq)
             
             ## BREAD
             # diagonal
             breadInv <- as.matrix(Matrix::bdiag(lapply(fit, function(x) solve(sandwich::bread(x) * nrow(x$model) / sum(summary(x)$df[1:2])))))
             dimnames(breadInv) <- list(dimnames, dimnames)
             
             # off-diagonal        
             # deriv12
             if (inherits(expData, "weightData")) {
               derivFUN <- switch(fit[[2]]$family$family, 
                                  gaussian = deriv(~ (- (M - mu)^2 / (2 * sigma^2)), "mu"), 
                                  binomial = deriv(~ (M * log(mu) + (1-M) * log(1-mu)), "mu"), 
                                  poisson = deriv(~ (M * log(mu) - mu), "mu"))
               sigma <- sqrt(summary(fit[[2]])$dispersion)
               # first
               X12 <- adaptx(expData, fit[[1]], obs = FALSE)
               mu <- predict(fit[[2]], newdat = X12$newdat, type = "response")
               M <- X12$newdat[, attr(terms(expData), "vartype")$M]
               if (!is.numeric(M)) M <- as.numeric(M) - 1
               deriv12a <- X12$modmat * fit[[2]]$family$mu.eta(predict(fit[[2]], newdat = X12$newdat)) * as.vector(attr(eval(derivFUN), "gradient"))
               
               # second
               X12 <- adaptx(expData, fit[[1]], obs = TRUE)
               mu <- predict(fit[[2]], newdat = X12$newdat, type = "response")
               M <- X12$newdat[, attr(terms(expData), "vartype")$M]
               if (!is.numeric(M)) M <- as.numeric(M) - 1
               deriv12b <- X12$modmat * fit[[2]]$family$mu.eta(predict(fit[[2]], newdat = X12$newdat)) * as.vector(attr(eval(derivFUN), "gradient"))
               
               deriv12 <- deriv12a - deriv12b
               breadInv[ind[[1]], ind[[2]]] <- -t(sandwich::estfun(fit[[1]])) %*% deriv12 / nrow(deriv12)
             }
             else if (inherits(expData, "impData")) {
               X12 <- adaptx(expData, fit[[1]], obs = FALSE)
               deriv12 <- X12$modmat * fit[[2]]$family$mu.eta(predict(fit[[2]], newdat = X12$newdat))
               
               breadInv[ind[[1]], ind[[2]]] <- -t(sandwich::estfun(fit[[1]]) / resid(fit[[1]], type = "response")) %*% deriv12 / nrow(deriv12)
             }
             
             # deriv13
             if (length(fit) > 2) {
               derivFUN <- switch(fit[[3]]$family$family, 
                                  gaussian = deriv(~ (- (X - mu)^2 / (2 * sigma^2)), "mu"), 
                                  binomial = deriv(~ (X * log(mu) + (1-X) * log(1-mu)), "mu"), 
                                  poisson = deriv(~ (X * log(mu) - mu), "mu"))
               sigma <- sqrt(summary(fit[[3]])$dispersion)
               X13 <- adaptx(expData, fit[[1]], obs = TRUE, xFit = fit[[3]])
               mu <- predict(fit[[3]], newdata = X13$newdat, type = "response")
               X <- X13$newdat[, attr(terms(expData), "vartype")$X]
               if (!is.numeric(X)) X <- as.numeric(X) - 1
               deriv13 <- - X13$modmat * fit[[3]]$family$mu.eta(predict(fit[[3]], newdat = X13$newdat)) * as.vector(attr(eval(derivFUN), "gradient"))
               breadInv[ind[[1]], ind[[3]]] <- -t(sandwich::estfun(fit[[1]])) %*% deriv13 / nrow(deriv13)
             }
    
             bread <- solve(breadInv)
             vcov <- as.matrix((bread %*% meat %*% t(bread)) / nrow(estEq))
             
             neModel <- list(neModelFit = neModelFit, vcov = vcov[ind[[1]], ind[[1]]])
           })
    terms <- neModelFit$terms
    attr(terms, "vartype") <- attr(attr(expData, "terms"), "vartype")
    class(terms) <- c("neTerms.object", class(terms))
    neModel <- c(neModel, terms = terms)
    class(neModel) <- c(if (se == "bootstrap") "neModelBoot", "neModel")
    return(neModel)
}

#' @rdname plot.neModel
#' @export
plot.neModel <- function (x, xRef, covLev, level = 0.95, 
    transf = identity, ylabels, yticks.at, ...) 
{
    args <- as.list(match.call())
    args[[1]] <- substitute(plot)
    args$x <- substitute(neEffdecomp(x, xRef = xRef, covLev = covLev))
    effdecomp <- eval(args$x)
    args[c("xRef", "covLev")] <- NULL
    eval(as.call(args))
    effdecompMessage <- if (substitute(transf) == "identity") "Effect decomposition on the scale of the linear predictor\n" else paste0("Effect decomposition on ", substitute(transf), "(scale of the linear predictor)\n")
    covModifier <- dimnames(attr(effdecomp, "covLev"))[[2]] %in% attr(effdecomp, "covModifier")
    sep <- if (all(covModifier) | !any(covModifier)) "" else ", "
    covMessage <- if (!is.null(attr(effdecomp, "covLev"))) paste("conditional on:", paste(paste(dimnames(attr(effdecomp, "covLev")[, covModifier, drop = FALSE])[[2]], attr(effdecomp, "covLev")[, covModifier, drop = FALSE], sep = " = ", collapse = ", "), 
                                                                                          paste(dimnames(attr(effdecomp, "covLev")[, !covModifier, drop = FALSE])[[2]], collapse = ", "), sep = sep), "\n") 
    message(effdecompMessage, covMessage,
            paste0("with x* = ", attr(effdecomp, "xRef")[1], ", x = ", attr(effdecomp, "xRef")[2]), "\n")   
}

#' @rdname plot.neModel
#' @export
plot.neModelBoot <- function (x, xRef, covLev, level = 0.95, 
    ci.type = "norm", transf = identity, ylabels, yticks.at, ...) 
{
    args <- as.list(match.call())
    args[[1]] <- substitute(plot)
    args$x <- substitute(neEffdecomp(x, xRef = xRef, covLev = covLev))
    effdecomp <- eval(args$x)
    args[c("xRef", "covLev")] <- NULL
    eval(as.call(args))
    effdecompMessage <- if (substitute(transf) == "identity") "Effect decomposition on the scale of the linear predictor\n" else paste0("Effect decomposition on ", substitute(transf), "(scale of the linear predictor)\n")
    covModifier <- dimnames(attr(effdecomp, "covLev"))[[2]] %in% attr(effdecomp, "covModifier")
    sep <- if (all(covModifier) | !any(covModifier)) "" else ", "
    covMessage <- if (!is.null(attr(effdecomp, "covLev"))) paste("conditional on:", paste(paste(dimnames(attr(effdecomp, "covLev")[, covModifier, drop = FALSE])[[2]], attr(effdecomp, "covLev")[, covModifier, drop = FALSE], sep = " = ", collapse = ", "), 
                                                                                          paste(dimnames(attr(effdecomp, "covLev")[, !covModifier, drop = FALSE])[[2]], collapse = ", "), sep = sep), "\n") 
    message(effdecompMessage, covMessage,
            paste0("with x* = ", attr(effdecomp, "xRef")[1], ", x = ", attr(effdecomp, "xRef")[2]), "\n") 
}

#' @method print neBootCI
#' @export
print.neBootCI <- function (x, ...) 
{
    attr <- attributes(x)
    attributes(x)[c("level", "coef", "R", "type", "class")] <- NULL
    print.default(x)
    cat("---\n")
    cat(switch(attr$type, norm = "normal approximation bootstrap", 
        basic = "basic bootstrap", stud = "studentized bootstrap", 
        perc = "bootstrap percentile", bca = "adjusted bootstrap percentile (BCa)"), 
        " confidence intervals\nbased on ", attr$R, " bootstrap samples", 
        sep = "")
    cat("\n\n")
}

#' @method print neModelCI
#' @export
print.neModelCI <- function (x, ...) 
{
    attributes(x)[c("level", "coef", "class")] <- NULL
    print.default(x)
}

#' @method print summary.neModel
#' @export
print.summary.neModel <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("Natural effect model\n")
    if (attr(x, "class_object")[1] == "neModelBoot") cat("with standard errors based on the non-parametric bootstrap\n---\n")
    else if (attr(x, "class_object")[1] == "neModel") cat("with robust standard errors based on the sandwich estimator\n---\n")
    cat("Exposure:", x$terms$X, "\nMediator(s):", paste(x$terms$M, 
        collapse = ", "), "\n---\n")
    cat("Parameter estimates:\n")
    printCoefmat(x$coef.table, digits = digits, has.Pvalue = TRUE, 
        P.values = TRUE)
}

#' @rdname expData-methods
#' @export
residuals.expData <- function(object, ...) 
{
    residuals(attr(object, "model"), ...)
}

#' @rdname expData-methods
#' @export
residualPlot.expData <- function(object, ...) 
{
    car::residualPlot(attr(object, "model"), ...)
}

#' @rdname expData-methods
#' @export
residualPlots.expData <- function(object, ...) 
{
    car::residualPlots(attr(object, "model"), ...)
}

#' @rdname neModel-methods
#' @export
residualPlot.neModel <- function(object, ...) 
{
    car::residualPlot(object$neModelFit, ...)
}

#' @rdname neModel-methods
#' @export
residualPlots.neModel <- function(object, ...) 
{
    object$neModelFit$call[[1]] <- quote(glm)
    object$neModelFit$call[["data"]] <- object$neModelFit$call[["expData"]]
    object$neModelFit$call[["expData"]] <- NULL
    object$neModelFit$call[["weights"]] <- quote(weights(object$neModelFit, type = "prior"))
    if (!is.null(object$neModelFit$call[["xFit"]])) object$neModelFit$call[["xFit"]] <- NULL
    object$neModelFit <- update(object$neModelFit)
    car::residualPlots(object$neModelFit, ...)
}

#' @rdname neModel-methods
#' @method summary neModel
#' @export
summary.neModel <- function (object, ...) 
{
    summary <- summary(object$neModelFit)
    se <- sqrt(diag(vcov(object)))
    zvalue <- coef(object)/se
    pvalue <- 2 * pnorm(-abs(zvalue))
    coef.table <- cbind(coef(object), se, zvalue, pvalue)
    dimnames(coef.table) <- list(names(coef(object)), c("Estimate", 
        "Std. Error", "z value", "Pr(>|z|)"))
    summary$coef.table <- coef.table
    summary$terms <- attr(object$terms, "vartype")
    class(summary) <- "summary.neModel"
    attr(summary, "class_object") <- class(object) 
    return(summary)
}

#' @rdname neModel-methods
#' @export
vcov.neModel <- function (object, ...) 
{
    return(object$vcov)
}

#' @export
vcov.neModelBoot <- function (object, ...) 
{
    covmat <- cov(object$bootRes$t)
    dimnames(covmat) <- dimnames(summary.glm(object$neModelFit)$cov.scaled)
    return(covmat)
}

#' @rdname expData-methods
#' @export
weights.expData <- function (object, ...)
{
    attr(object, "weights")
}

#' @rdname neModel-methods
#' @export
weights.neModel <- function (object, ...) 
{
    weights(object$neModelFit, ...)
}