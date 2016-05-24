#' Expand the dataset and calculate ratio-of-mediator probability weights
#'
#' @description This function both expands the data along hypothetical exposure values and calculates ratio-of-mediator probability weights.
#' @param object an object used to select a method.
#' @param ... additional arguments.
#' @return A data frame of class \code{c("data.frame", "expData", "weightData")}. See \code{\link{expData}} for its structure.
#' @details Generic function that both expands the data along hypothetical exposure values and
#' calculates ratio-of-mediator probability weights 
#' 
#' \deqn{\frac{\hat P(M_i \vert X_i = x^*, C_i)}{\hat P(M_i \vert X_i = x, C_i)}}
#' 
#' for each observation unit \emph{i} in this expanded dataset in a single run.
#' These weights are ratios of probabilities or probability densities from the mediator model distribution, which can be specified either externally as a fitted model object (\code{\link{neWeight.default}})
#' or internally (\code{\link{neWeight.formula}}).
#' @seealso \code{\link{neWeight.default}}, \code{\link{neWeight.formula}}, \code{\link{expData}}
#' @references
#' Hong, G. (2010). Ratio of mediator probability weighting for estimating natural direct and indirect effects. In \emph{Proceedings of the American Statistical Association, Biometrics Section}, pp. 2401-2415. American Statistical Association, Alexandria, VA.
#' 
#' Lange, T., Vansteelandt, S., & Bekaert, M. (2012). A Simple Unified Approach for Estimating Natural Direct and Indirect Effects. \emph{American Journal of Epidemiology}, \bold{176}(3), 190-195.
#' @export
neWeight <- function (object, ...) 
{
    UseMethod("neWeight")
}

#' Expand the dataset and calculate ratio-of-mediator probability weights
#'
#' @description This function both expands the data along hypothetical exposure values and calculates ratio-of-mediator probability weights.
#' @param object fitted model object representing the mediator model.
#' @param formula a \code{\link[stats]{formula}} object providing a symbolic description of the mediator model. Redundant if already specified in call for fitted model specified in \code{object} (see details).
#' @inheritParams neImpute.default
#' @return A data frame of class \code{c("data.frame", "expData", "weightData")}. See \code{\link{expData}} for its structure.
#' @details The calculated weights are ratios of fitted probabilities or probability densities from the distribution of the mediator model.
#' This model needs to be specified as a fitted object in the \code{object} argument.
#'
#' If the model-fitting function used to fit the mediator model does not require specification of a \code{formula} or \code{data} argument,
#' these need to be specified explicitly in order to enable \code{neWeight.default} to extract pointers to variable types relevant for mediation analysis.
# (also see \code{\link{neTerms}}).
#'
#' Whether a \code{\link[stats]{formula}} is specified externally (in the call for the fitted mediator model object which is specified in \code{object}) or internally (via the \code{formula} argument),
#' it always needs to be of the form \code{M ~ X + C1 + C2}, with predictor variables entered in the following prespecified order:
#' \enumerate{
#'  \item exposure \code{X}: The first predictor is coded as exposure or treatment.
#'  \item baseline covariates \code{C}: All remaining predictor variables are automatically coded as baseline covariates.
#' }
#'
#' It is important to adhere to this prespecified order to enable \code{neWeight} to create valid pointers to these different types of predictor variables.
#' This requirement extends to the use of operators different than the \code{+} operator, such as the \code{:} and \code{*} operators (when e.g. adding interaction terms). 
#' For instance, the formula specifications \code{M ~ X * C1 + C2}, \code{M ~ X + C1 + X:C1 + C2} and \code{Y ~ X + X:C1 + C1 + C2} will create identical pointers to the different types of variables,
#' as the order of the unique predictor variables is identical in all three specifications. 
#' 
#' Furthermore, categorical exposures that are not coded as factors in the original dataset, should be specified as factors in the formula, 
#' using the \code{\link[base]{factor}} function, e.g. \code{M ~ factor(X) + C1 + C2}. 
#' Quadratic or higher-order polynomial terms can be included as well, by making use of the \code{\link[base]{I}} function or by using the \code{\link[stats]{poly}} function.
#' For instance, \code{M ~ X + I(X^2) + C1 + C2} and \code{M ~ poly(X, 2, raw = TRUE) + C1 + C2} are equivalent and result in identical pointers to the different types of variables.
# We do not recommend the use of orthogonal polynomials (i.e. using the default argument specification \code{raw = FALSE} in \code{poly}).
#' 
#' The command \code{terms(object, "vartype")} (with \code{object} replaced by the name of the resulting expanded dataset) can be used to check whether valid pointers have been created.
#'
#' In contrast to imputation models with categorical exposures, additional arguments need to be specified if the exposure is continuous.
#' All of these additional arguments are related to the sampling procedure for the exposure.
#'
#' Whereas the number of replications \code{nRep} for categorical variables equals the number of levels for the exposure coded as a factor (i.e. the number of hypothetical exposure values), the number of desired replications needs to be specified explicitly for continuous exposures.
#' Its default is 5.
#'
#' If \code{xFit} is left unspecified, the hypothetical exposure levels are automatically sampled from a linear model for the exposure, conditional on a linear combination of all covariates.
#' If one wishes to use another model for the exposure, this default model specification can be overruled by referring to a fitted model object in the \code{xFit} argument.
#' Misspecification of this sampling model does not induce bias in the estimated coefficients and standard errors of the natural effect model.
#'
#' The \code{xSampling} argument allows to specify how the hypothetical exposure levels should be sampled from the conditional exposure distribution (which is either entered explicitly using the \code{xFit} argument or fitted automatically as described in the previous paragraph).
#' The \code{"random"} option randomly samples \code{nRep} draws from the exposure distribution, whereas the \code{"quantiles"} option (default) samples \code{nRep} quantiles at equal-sized probability intervals. Only the latter hence yields fixed exposure levels given \code{nRep} and \code{xFit}. \cr\cr
#' In order to guarantee that the entire support of the distribution is being sampled (which might be a concern if \code{nRep} is chosen to be small), the default lower and upper sampled quantiles are the 5th and 95th percentiles.
#' The intermittent quantiles correspond to equal-sized probability intervals. So, for instance, if \code{nRep = 4}, then the sampled quantiles will correspond to probabilities 0.05, 0.35, 0.65 and 0.95.
#' These default 'outer' quantiles can be changed by specifying the \code{percLim} argument accordingly. By specifying \code{percLim = NULL}, the standard quantiles will be sampled (e.g., 0.2, 0.4, 0.6 and 0.8 if \code{nRep = 4}).
#'
#' @seealso \code{\link{neWeight}}, \code{\link{neWeight.formula}}, \code{\link{expData}}
#' @examples
#' data(UPBdata)
#' 
#' ## example using glm
#' fit.glm <- glm(negaff ~ att + gender + educ + age, data = UPBdata)
#' weightData <- neWeight(fit.glm, nRep = 2)
#' 
#' ## example using vglm (yielding identical results as with glm)
#' library(VGAM)
#' fit.vglm <- vglm(negaff ~ att + gender + educ + age, 
#'                  family = gaussianff, data = UPBdata)
#' weightData2 <- neWeight(fit.vglm, nRep = 2)
#'
#' \dontshow{
#' library(VGAM) 
#' expData <- neWeight(negaff ~ factor(attbin) + gender + educ + age, family = gaussianff, data = UPBdata, FUN = vglm)
#' neMod <- neModel(UPB ~ attbin0 + attbin1 + gender + educ + age, family = binomial, expData = expData, nBoot = 2)
#'  
#' fit1 <- glm(negaff ~ att + gender + educ + age, data = UPBdata)
#' expData1 <- neWeight(fit1)
#' w1 <- attr(expData1, "weights")
#' expData1f <- neWeight(negaff ~ att + gender + educ + age, data = UPBdata)
#' w1f <- attr(expData1f, "weights")
#' head(expData1); head(expData1f)
#' head(w1); head(w1f)
#'
#' # test vglm (vglm is also vgam class, but not other way around!)
#' fit1b <- vgam(negaff ~ att + gender + educ + age, family = gaussianff, data = UPBdata)
#' expData1b <- neWeight(fit1b)
#' head(attr(expData1, "weights")); head(attr(expData1b, "weights"))
#' fit1b <- vgam(negaff ~ s(att) + gender + educ + age, family = gaussianff, data = UPBdata)
#' expData1b <- neWeight(fit1b)
#' head(attr(expData1, "weights")); head(attr(expData1b, "weights"))
#' expData1bf <- neWeight(negaff ~ s(att) + gender + educ + age, FUN = vgam, family = gaussianff, data = UPBdata)
#' head(attr(expData1b, "weights")); head(attr(expData1bf, "weights"))
#' ##
#'
#' UPBdata$negaff2 <- cut(UPBdata$negaff, breaks = 2, labels = c("low", "high"))
#' fit2 <- glm(negaff2 ~ att + gender + educ + age, family = binomial, data = UPBdata)
#' expData2 <- neWeight(fit2)
#' w2 <- attr(expData2, "weights")
#' expData2f <- neWeight(negaff2 ~ att + gender + educ + age, family = binomial, data = UPBdata)
#' w2f <- attr(expData2f, "weights")
#' head(expData2); head(expData2f)
#' head(w2); head(w2f)
#'
#' # test vglm
#' fit2b <- vgam(negaff2 ~ att + gender + educ + age, family = binomialff, data = UPBdata)
#' expData2b <- neWeight(fit2b)
#' head(attr(expData2, "weights")); head(attr(expData2b, "weights"))
#' fit2b <- vgam(negaff2 ~ s(att) + gender + educ + age, family = binomialff, data = UPBdata)
#' expData2b <- neWeight(fit2b)
#' head(attr(expData2, "weights")); head(attr(expData2b, "weights"))
#' expData2bf <- neWeight(negaff2 ~ s(att) + gender + educ + age, FUN = vgam, family = binomialff, data = UPBdata)
#' head(attr(expData2b, "weights")); head(attr(expData2bf, "weights"))
#' ##
#'
#' UPBdata$negaff3 <- cut(UPBdata$negaff, breaks = 3, labels = c("low", "moderate", "high"))
#' UPBdata$negaff3 <- as.numeric(UPBdata$negaff3)
#' fit3 <- glm(negaff3 ~ att + gender + educ + age, family = "poisson", data = UPBdata)
#' expData3 <- neWeight(fit3)
#' w3 <- attr(expData3, "weights")
#' expData3f <- neWeight(negaff3 ~ att + gender + educ + age, family = poisson, data = UPBdata)
#' w3f <- attr(expData3f, "weights")
#' head(expData3); head(expData3f)
#' head(w3); head(w3f)
#'
#' # test vglm
#' fit3b <- vgam(negaff3 ~ att + gender + educ + age, family = poissonff, data = UPBdata)
#' expData3b <- neWeight(fit3b)
#' head(attr(expData3, "weights")); head(attr(expData3b, "weights"))
#' fit3b <- vgam(negaff3 ~ s(att) + gender + educ + age, family = poissonff, data = UPBdata)
#' expData3b <- neWeight(fit3b)
#' head(attr(expData3, "weights")); head(attr(expData3b, "weights"))
#' expData3bf <- neWeight(negaff3 ~ s(att) + gender + educ + age, FUN = vgam, family = poissonff, data = UPBdata)
#' head(attr(expData3b, "weights")); head(attr(expData3bf, "weights"))}
#' @export
neWeight.default <- function (object, formula, data, nRep = 5, xSampling = c("quantiles", 
    "random"), xFit, percLim = c(0.05, 0.95), ...) 
{
    args <- as.list(match.call())[-1L]
    if (!is.null(args$nMed) && args$nMed != 1) warning("The joint mediation approach is only available for the imputation-based approach! nMed = 1 was specified instead.")
    nMed <- args$nMed <- 1
    fit <- object
    args$data <- if (missing(data)) {
        extrCall(fit)$data
    }
    else substitute(data)
    if (!isTRUE(args$skipExpand)) {
        if (missing(formula)) 
            formula <- extrCall(fit)$formula
        class(formula) <- c("Mformula", "formula")
        vartype <- args$vartype <- attr(neTerms(formula, Y = NA, 
            nMed = nMed, joint = joint), "vartype")
        xFact <- if ("SuperLearner" %in% class(fit)) {
            is.factor(model.frame(formula, eval(args$data))[, attr(vartype, "xasis")])
        } else if (any(c("vglm", "vgam") %in% class(fit))) {
            is.factor(VGAM::model.frame(fit)[, attr(vartype, "xasis")])
        } else {
            is.factor(model.frame(fit)[, attr(vartype, "xasis")])
        }
        if (xFact) {
            dontapply <- c("nRep", "xSampling", "xFit", "percLim")
            ind <- !sapply(args[dontapply], is.null)
            ind <- ind[which(ind)]
            if (length(ind) > 1) 
                warning(gettextf("The arguments %s don't apply when the exposure is categorical!", 
                  paste(paste0("'", names(ind), "'"), collapse = ", ")))
            else if (length(ind) == 1) 
                warning(gettextf("The argument '%s' doesn't apply when the exposure is categorical!", 
                  names(ind)))
        }
        args$xSampling <- xSampling <- match.arg(xSampling)
        if (all(!missing(percLim), xSampling == "random")) 
            warning("The percentile limits that have been specified in 'percLim' do not apply when 'xSampling' is specified as 'random'!")
        joint <- args$joint <- TRUE
        if (missing(nRep)) 
            args$nRep <- substitute(nRep)
        if (missing(percLim)) 
            args$percLim <- percLim
        expData <- do.call("expandData", c(x = substitute(vartype$X), 
            args))
        nExp <- ifelse(joint, 1, nMed)
        tmp <- vartype$Xexp <- paste0(vartype$X, seq(0, nExp))
        names(expData)[sapply(paste0("aux", seq(0, nExp)), grep, 
            names(expData))] <- tmp
        other <- names(expData)[!is.element(names(expData), c("id", 
            vartype$X, vartype$Xexp))]
        expData <- expData[, c("id", vartype$Xexp, other)]
    }
    else {
        if (is.null(args$expData)) 
            stop("expData is missing!")
        expData <- eval.parent(args$expData)
        attr <- attributes(expData)
        vartype <- attr(attr$terms, "vartype")
    }
    family <- if(inherits(fit, "vglm")) fit@family@vfamily[1] else fit$family$family
    family <- c("gaussian", "binomial", "poisson", "multinomial")[mapply(function(x, 
        y) grepl(y, x), as.character(family), c("gaussian", "binomial", 
        "poisson", "multinomial"))]
    dispersion <- if (inherits(fit, "vglm")) 
        fit@misc$dispersion
    else summary(fit)$dispersion
    dfun <- switch(family, 
         gaussian = function(x) dnorm(x, mean = predict(fit, newdata = expData, type = "response"), sd = sqrt(dispersion)), 
         binomial = function(x) {if (is.factor(x)) x <- as.numeric(x) - 1
                     return(dbinom(x, size = 1, prob = predict(fit, newdata = expData, type = "response")))}, 
         poisson = function(x) dpois(x, lambda = predict(fit, newdata = expData, type = "response")),
         multinomial = function(x) {pred <- predict(fit, newdata = expData, type = "response")
                        return(sapply(1:nrow(expData), function(i) pred[i, as.character(x[i])]))})
    expData[, vartype$X] <- expData[, vartype$Xexp[2]]
    try(weightsNum <- dfun(expData[, vartype$M]), silent = TRUE)
    checkExist <- exists("weightsNum")
    if (!checkExist) expData[, attr(vartype, "xasis")] <- expData[, vartype$Xexp[2]] 
    weightsNum <- dfun(expData[, vartype$M])
    expData[, vartype$X] <- expData[, vartype$Xexp[1]]
    if (!checkExist) expData[, attr(vartype, "xasis")] <- expData[, vartype$Xexp[1]] 
    weightsDenom <- dfun(expData[, vartype$M])
    expData <- expData[, -ncol(expData)]
    if (!isTRUE(args$skipExpand)) {
        attributes(expData) <- c(attributes(expData), list(model = fit, 
            call = match.call(), terms = neTerms(formula, nMed, 
                joint), weights = weightsNum/weightsDenom))
        attr(attr(expData, "terms"), "vartype") <- vartype
        class(expData) <- c(class(expData), "expData", "weightData")
    }
    else {
        attributes(expData) <- attr
        attr(expData, "model") <- fit
        attr(expData, "weights") <- weightsNum/weightsDenom
    }
    return(expData)
}

#' Expand the dataset and calculate ratio-of-mediator probability weights
#'
#' @description This function both expands the data along hypothetical exposure values and calculates ratio-of-mediator probability weights.
#' @param object a \code{\link[stats]{formula}} object providing a symbolic description of the mediator model (see details).
#' @inheritParams neImpute.formula
#' @return A data frame of class \code{c("data.frame", "expData", "weightData"))}. See \code{\link{expData}} for its structure.
#' @details The calculated weights are ratios of fitted probabilities or probability densities from the distribution of the mediator model.
#' This model is fitted internally by extracting information from the arguments \code{object}, \code{family}, \code{data}, \code{FUN} and \code{...}.
#'
#' For mediation model specification via the \code{object} argument, use a \code{\link[stats]{formula}} of the form\cr\code{M ~ X + C1 + C2},
#' with predictor variables entered in the following prespecified order:
#' \enumerate{
#'  \item exposure \code{X}: The first predictor is coded as exposure or treatment.
#'  \item baseline covariates \code{C}: All remaining predictor variables are automatically coded as baseline covariates.
#' }
#'
#' It is important to adhere to this prespecified order to enable \code{neWeight} to create valid pointers to these different types of predictor variables.
#' This requirement extends to the use of operators different than the \code{+} operator, such as the \code{:} and \code{*} operators (when e.g. adding interaction terms). 
#' For instance, the formula specifications \code{M ~ X * C1 + C2}, \code{M ~ X + C1 + X:C1 + C2} and \code{Y ~ X + X:C1 + C1 + C2} will create identical pointers to the different types of variables,
#' as the order of the unique predictor variables is identical in all three specifications. 
#' 
#' Furthermore, categorical exposures that are not coded as factors in the original dataset, should be specified as factors in the formula, 
#' using the \code{\link[base]{factor}} function, e.g. \code{M ~ factor(X) + C1 + C2}. 
#' Quadratic or higher-order polynomial terms can be included as well, by making use of the \code{\link[base]{I}} function or by using the \code{\link[stats]{poly}} function.
#' For instance, \code{M ~ X + I(X^2) + C1 + C2} and \code{M ~ poly(X, 2, raw = TRUE) + C1 + C2} are equivalent and result in identical pointers to the different types of variables.
# We do not recommend the use of orthogonal polynomials (i.e. using the default argument specification \code{raw = FALSE} in \code{poly}).
#' 
#' The command \code{terms(object, "vartype")} (with \code{object} replaced by the name of the resulting expanded dataset) can be used to check whether valid pointers have been created.
#'
#' The type of mediator model can be defined by specifying an appropriate model-fitting function via the \code{FUN} argument (its default is \code{\link[stats]{glm}}).
#' This method can only be used with model-fitting functions that require a \code{formula} argument.
#'
#' In contrast to imputation models with categorical exposures, additional arguments need to be specified if the exposure is continuous.
#' All of these additional arguments are related to the sampling procedure for the exposure.
#'
#' Whereas the number of replications \code{nRep} for categorical variables equals the number of levels for the exposure coded as a factor (i.e. the number of hypothetical exposure values), the number of desired replications needs to be specified explicitly for continuous exposures.
#' Its default is 5.
#'
#' If \code{xFit} is left unspecified, the hypothetical exposure levels are automatically sampled from a linear model for the exposure, conditional on a linear combination of all covariates.
#' If one wishes to use another model for the exposure, this default model specification can be overruled by referring to a fitted model object in the \code{xFit} argument.
#' Misspecification of this sampling model does not induce bias in the estimated coefficients and standard errors of the natural effect model.
#'
#' The \code{xSampling} argument allows to specify how the hypothetical exposure levels should be sampled from the conditional exposure distribution (which is either entered explicitly using the \code{xFit} argument or fitted automatically as described in the previous paragraph).
#' The \code{"random"} option randomly samples \code{nRep} draws from the exposure distribution, whereas the \code{"quantiles"} option (default) samples \code{nRep} quantiles at equal-sized probability intervals. Only the latter hence yields fixed exposure levels given \code{nRep} and \code{xFit}. \cr\cr
#' In order to guarantee that the entire support of the distribution is being sampled (which might be a concern if \code{nRep} is chosen to be small), the default lower and upper sampled quantiles are the 5th and 95th percentiles.
#' The intermittent quantiles correspond to equal-sized probability intervals. So, for instance, if \code{nRep = 4}, then the sampled quantiles will correspond to probabilities 0.05, 0.35, 0.65 and 0.95.
#' These default 'outer' quantiles can be changed by specifying the \code{percLim} argument accordingly. By specifying \code{percLim = NULL}, the standard quantiles will be sampled (e.g., 0.2, 0.4, 0.6 and 0.8 if \code{nRep = 4}).
#'
#' @seealso \code{\link{neWeight.default}}, \code{\link{expData}}
#' @examples
#' data(UPBdata)
#' 
#' ## example using glm
#' weightData <- neWeight(negaff ~ att + gender + educ + age, 
#'                        data = UPBdata, nRep = 2)
#' 
#' ## example using vglm (yielding identical results as with glm)
#' library(VGAM)
#' weightData2 <- neWeight(negaff ~ att + gender + educ + age, 
#'                         family = gaussianff, data = UPBdata, nRep = 2, FUN = vglm)
#' @export
neWeight.formula <- function (object, family, data, FUN = glm, nRep = 5, xSampling = c("quantiles", 
    "random"), xFit, percLim = c(0.05, 0.95), ...) 
{
    args <- as.list(match.call())[-1L]
    if (missing(data)) 
      data <- environment(object)
    formula <- object
    if (missing(family)) 
      family <- formals(FUN)$family
    argsFUN <- list(formula = formula, family = family, data = data, ...)
    args$object <- do.call(FUN, eval(argsFUN[!names(argsFUN) %in% "nMed"]))  
    call <- substitute(list(formula = formula, family = family, data = data, ...))
    call[[1]] <- substitute(FUN)
    call$nMed <- NULL
    if (isS4(args$object)) 
        args$object@call <- call
    else args$object$call <- call
    args$formula <- args$family <- args$data <- NULL
    expData <- do.call("neWeight.default", args)
    return(expData)
}