#' Expand the dataset and impute nested counterfactual outcomes
#'
#' @description This function both expands the data along hypothetical exposure values and imputes nested counterfactual outcomes.
#' @param object an object used to select a method.
#' @param ... additional arguments.
#' @return A data frame of class \code{c("data.frame", "expData", "impData")}. See \code{\link{expData}} for its structure.
#' @details Generic function that both expands the data along hypothetical exposure values (for each observation unit \emph{i}) and imputes nested counterfactual outcomes in this expanded dataset in a single run.
#' Imputed counterfactual outcomes 
#' 
#' \deqn{\hat E(Y_i \vert X_i = x, M_i, C_i)}
#' 
#' are predictions from the imputation model that can be specified either externally as a fitted model object (\code{\link{neImpute.default}})
#' or internally (\code{\link{neImpute.formula}}).
#' @seealso \code{\link{neImpute.default}}, \code{\link{neImpute.formula}}, \code{\link{neModel}}, \code{\link{expData}}
#' @references
#' Vansteelandt, S., Bekaert, M., & Lange, T. (2012). Imputation Strategies for the Estimation of Natural Direct and Indirect Effects. \emph{Epidemiologic Methods}, \bold{1}(1), Article 7.
#' 
#' Loeys, T., Moerkerke, B., De Smet, O., Buysse, A., Steen, J., & Vansteelandt, S. (2013). Flexible Mediation Analysis in the Presence of Nonlinear Relations: Beyond the Mediation Formula. \emph{Multivariate Behavioral Research}, \bold{48}(6), 871-894.
#' @export
neImpute <- function (object, ...) 
{
    UseMethod("neImpute")
}

#' Expand the dataset and impute nested counterfactual outcomes
#'
#' @description This function both expands the data along hypothetical exposure values and imputes nested counterfactual outcomes.
#' @param object fitted model object representing the imputation model.
#' @param formula a \code{\link[stats]{formula}} object providing a symbolic description of the imputation model. Redundant if already specified in call for fitted model specified in \code{object} (see details).
#' @param data data, as matrix or data frame, containing the exposure (and other relevant) variables. Redundant if already specified in call for fitted model specified in \code{object} (see details).
#' @param nMed number of mediators.
#' @param nRep number of replications or hypothetical values of the exposure to sample for each observation unit.
#' @param xSampling character string indicating how to sample from the conditional exposure distribution.
#' Possible values are \code{"quantiles"} or \code{"random"} (see details).
#' @param xFit an optional fitted object (preferably \code{glm}) for the conditional exposure distribution (see details).
#' @param percLim a numerical vector of the form \code{c(lower, upper)} indicating the extreme percentiles to sample when using \code{"quantiles"} as sampling method to sample from the conditional exposure distribution (see details).
#' @param ... additional arguments.
#' @return A data frame of class \code{c("data.frame", "expData", "impData")}. See \code{\link{expData}} for its structure.
#' @details Imputed counterfactual outcomes are predictions from the imputation model that needs to be specified as a fitted object in the \code{object} argument.
#'
#' If the model-fitting function used to fit the imputation model does not require specification of a \code{formula} or \code{data} argument (when using e.g. \code{\link[SuperLearner]{SuperLearner}}),
#' these need to be specified explicitly in order to enable \code{neImpute.default} to extract pointers to variable types relevant for mediation analysis.
#'
#' Whether a \code{\link[stats]{formula}} is specified externally (in the call for the fitted imputation model object which is specified in \code{object}) or internally (via the \code{formula} argument),
#' it always needs to be of the form \code{Y ~ X + M1 + M2 + M3 + C1 + C2}, with the same outcome as in the final natural effect model and with predictor variables entered in the following prespecified order:
#' \enumerate{
#'  \item exposure \code{X}: The first predictor is coded as exposure or treatment.
#'  \item mediator(s) \code{M}: The second predictor is coded as mediator. In case of multiple mediators (\code{nMed > 1}), then predictors \code{2:(nMed + 1)} are coded as mediators.
#: if one wishes to estimate path-specific effects wrt multiple mediators (including exposure-induced confounders) in the natural effect model (\code{joint = FALSE}), one has to make sure to enter the mediators in the appropriate chronological order.
#'  \item baseline covariates \code{C}: All remaining predictor variables are automatically coded as baseline covariates.
#' }
#' 
#' It is important to adhere to this prespecified order to enable \code{neImpute} to create valid pointers to these different types of predictor variables.
#' This requirement extends to the use of operators different than the \code{+} operator, such as the \code{:} and \code{*} operators (when e.g. adding interaction terms). 
#' For instance, the formula specifications \code{Y ~ X * M + C1 + C2}, \code{Y ~ X + M + X:M + C1 + C2} and \code{Y ~ X + X:M + M + C1 + C2} will create identical pointers to the different types of variables,
#' as the order of the unique predictor variables is identical in all three specifications. 
#' 
#' Furthermore, categorical exposures that are not coded as factors in the original dataset, should be specified as factors in the formula, 
#' using the \code{\link[base]{factor}} function, e.g. \code{Y ~ factor(X) + M + C1 + C2}. 
#' Quadratic or higher-order polynomial terms can be included as well, by making use of the \code{\link[base]{I}} function or by using the \code{\link[stats]{poly}} function.
#' For instance, \code{Y ~ X + I(X^2) + M + C1 + C2} and \code{Y ~ poly(X, 2, raw = TRUE) + M + C1 + C2} are equivalent and result in identical pointers to the different types of variables.
# We do not recommend the use of orthogonal polynomials (i.e. using the default argument specification \code{raw = FALSE} in \code{poly}).
#' 
#' The command \code{terms(object, "vartype")} (with \code{object} replaced by the name of the resulting expanded dataset) can be used to check whether valid pointers have been created.
#'
#' If multiple mediators are specified (\code{nMed > 1}), the natural indirect effect parameter in the natural effect model captures the joint mediated effect. That is, the effect of the exposure on the outcome via these mediators considered jointly. 
#' The remaining effect of the exposure on the outcome (not mediated through the specified mediators) is then captured by the natural indirect effect parameter.
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
#' @seealso \code{\link{neImpute}}, \code{\link{neImpute.formula}}, \code{\link{neModel}}, \code{\link{expData}}
#' @examples
#' data(UPBdata)
#' 
#' ## example using glm imputation model with binary exposure
#' fit.glm <- glm(UPB ~ attbin + negaff + gender + educ + age, 
#'                family = binomial, data = UPBdata)
#' impData <- neImpute(fit.glm)
#' head(impData)
#' 
#' ## example using glm imputation model with continuous exposure
#' fit.glm <- glm(UPB ~ att + negaff + gender + educ + age, 
#'                family = binomial, data = UPBdata)
#' impData <- neImpute(fit.glm, nRep = 2)
#' head(impData)
#' 
#' ## example using vglm (yielding identical results as with glm)
#' library(VGAM)
#' fit.vglm <- vglm(UPB ~ att + negaff + gender + educ + age, 
#'                  family = binomialff, data = UPBdata)
#' impData2 <- neImpute(fit.vglm, nRep = 2)
#' head(impData2)
#' 
#' \donttest{## example using SuperLearner
#' library(Matrix)
#' library(SuperLearner)
#' SL.library <- c("SL.glm", "SL.glm.interaction", "SL.rpart",
#'                 "SL.step", "SL.stepAIC", "SL.step.interaction",
#'                 "SL.bayesglm", "SL.glmnet")
#' pred <- c("att", "negaff", "gender", "educ", "age")
#' fit.SL <- SuperLearner(Y = UPBdata$UPB, X = subset(UPBdata, select = pred),
#'                        SL.library = SL.library, family = binomial())
#' impSL <- neImpute(fit.SL, 
#'                   formula = UPB ~ att + negaff + gender + educ + age, 
#'                   data = UPBdata)
#' head(impSL)}
#' \dontshow{
#' #library(VGAM) 
#' #expData <- neImpute(UPB ~ factor(attbin) + negaff + gender + educ + age, family = binomialff, data = UPBdata, FUN = vglm)
#' #neMod <- neModel(UPB ~ attbin0 + attbin1 + gender + educ + age, family = binomial, expData = expData, nBoot = 2)
#' 
# UPBdata$att2 <- ifelse(UPBdata$attbin == "H", 1, 0)
#' UPBdata$att2 <- UPBdata$attbin
#' impData <- neImpute(UPB ~ factor(att2) * negaff + gender + educ + age, family = binomial, data = UPBdata)
#' impFit1 <- neModel(UPB ~ att20 * att21 + gender + educ + age, family = binomial, expData = impData, nBoot = 2)
#' impFit2 <- neModel(UPB ~ factor(att20) * factor(att21) + gender + educ + age, family = binomial, expData = impData, nBoot = 2)
#' summary(impFit1)
#' summary(impFit2)
#' 
#' head(neImpute(UPB ~ att + negaff + gender + educ + age, data = UPBdata, weights = rep(1, nrow(UPBdata))))
#' fit1 <- vglm(UPB ~ att + negaff + gender + educ + age, family = binomialff, data = UPBdata)
#' head(neImpute(fit1))
#' head(neImpute(UPB ~ att + negaff + gender + educ + age, family = binomialff, data = UPBdata, FUN = vglm, weights = rep(1, nrow(UPBdata))))
#' head(neImpute(UPB ~ att + negaff + gender + educ + age, family = binomial, data = UPBdata, weights = rep(1, nrow(UPBdata))))
#' UPBdata$att <- factor(cut(UPBdata$att, 2), labels = c("low", "high"))
#' fit2 <- glm(UPB ~ att * negaff * gender * educ * age, family = binomial, data = UPBdata)
#' # head(neImpute(fit2, nMed = 3, joint = FALSE))
#' # head(neImpute(UPB ~ att * negaff * gender * educ * age, data = UPBdata, family = binomial, nMed = 3, joint = FALSE))
#' head(neImpute(UPB ~ att + negaff + gender + educ + age, data = UPBdata, family = binomial))
#' head(neImpute(UPB ~ att + negaff + gender + educ + age, data = UPBdata))
#'
#' # # test with vglm!
#' library(VGAM)
#' impFit <- vglm(UPB ~ att + negaff + gender + educ + age, family = binomialff, data = UPBdata)
#' terms(impFit)
#' # debug(neImpute)
#' impData <- neImpute(impFit)
#' impFit2 <- glm(UPB ~ att + negaff + gender + educ + age, family = binomial, data = UPBdata)
#' impData2 <- neImpute(impFit2)
#' head(impData); head(impData2)
#' # check!!
#'
#' # test with gam?
#' library(gam)
#' impFit4 <- gam(UPB ~ att + negaff + gender + educ + age, family = binomial, data = UPBdata)
#' impData4 <- neImpute(impFit4)
#' head(impData2); head(impData4)
#' # check!}
#' @export
neImpute.default <- function (object, formula, data, nMed = 1, nRep = 5, xSampling = c("quantiles", 
    "random"), xFit, percLim = c(0.05, 0.95), ...) 
{
    args <- as.list(match.call())[-1L]
    fit <- object
    if (!isTRUE(args$skipExpand)) {
        args$data <- if (missing(data)) {
            extrCall(fit)$data
        }
        else substitute(data)
        if (missing(formula)) 
            formula <- extrCall(fit)$formula
        class(formula) <- c("Yformula", "formula")
        vartype <- args$vartype <- attr(neTerms(formula, nMed, 
            joint), "vartype")
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
        args$data <- eval.parent(args$data, n = 2)
        expData <- do.call("expandData", c(x = substitute(vartype$X), 
            args))
        nExp <- ifelse(joint, 1, nMed)
        tmp <- vartype$Xexp <- paste0(vartype$X, seq(0, nExp))
        names(expData)[sapply(paste0("aux", seq(0, nExp)), grep, 
            names(expData))] <- c(tmp[2:1], if (nExp > 1) tmp[(nExp + 
            1):3])
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
    expData[, vartype$X] <- expData[, vartype$Xexp[1]]    
    if ("SuperLearner" %in% class(fit)) {
        newdata <- expData[, names(eval(fit$call$X))]
        type <- NULL
        ind <- "pred"
    }
    else {
        newdata <- expData
        type <- "response"
        ind <- seq.int(nrow(newdata))
    }
    try(pred <- predict(fit, newdata = newdata, type = type)[ind], silent = TRUE)
    checkExist <- exists("pred")
    if (!checkExist) newdata[, attr(vartype, "xasis")] <- newdata[, vartype$Xexp[1]]  
    expData[, vartype$Y] <- predict(fit, newdata = newdata, type = type)[ind]
    expData <- expData[, -ncol(expData)]
    if (!isTRUE(args$skipExpand)) {
        attributes(expData) <- c(attributes(expData), list(model = fit, 
            call = match.call(), terms = neTerms(formula, nMed, 
                joint), weights = rep(1L, nrow(expData))))
        attr(attr(expData, "terms"), "vartype") <- vartype
        class(expData) <- c(class(expData), "expData", "impData")
    }
    else {
        attributes(expData) <- attr
        attr(expData, "model") <- fit
        attr(expData, "weights") <- rep(1L, nrow(expData))
    }
    return(expData)
}

#' Expand the dataset and impute nested counterfactual outcomes
#'
#' @description This function both expands the data along hypothetical exposure values and imputes nested counterfactual outcomes.
#' @param object a \code{\link[stats]{formula}} object providing a symbolic description of the imputation model (see details).
#' @param family a description of the error distribution and link function to be used in the model. Consult the help files of the model-fitting function specified in \code{FUN} for more details on appropriate argument specification.
#' @param FUN function used to fit model specified in \code{formula}.
#' @param ... additional arguments (passed to \code{FUN}).
#' @inheritParams neImpute.default
#' @return A data frame of class \code{c("data.frame", "expData", "impData")}. See \code{\link{expData}} for its structure.
#' @details Imputed counterfactual outcomes are predictions from the imputation model that is fitted internally by extracting information from the arguments \code{object}, \code{family}, \code{data}, \code{FUN} and \code{...}.
#'
#' For imputation model specification via the \code{object} argument, use a \code{\link[stats]{formula}} of the form 
#' 
#' \code{Y ~ X + M1 + M2 + M3 + C1 + C2},
#' 
#' with the same outcome as in the final natural effect model and with predictor variables entered in the following prespecified order:
#' \enumerate{
#'  \item exposure \code{X}: The first predictor is coded as exposure or treatment.
#'  \item mediator(s) \code{M}: The second predictor is coded as mediator. In case of multiple mediators (\code{nMed > 1}), then predictors \code{2:(nMed + 1)} are coded as mediators.
#'  \item baseline covariates \code{C}: All remaining predictor variables are automatically coded as baseline covariates.
#' }
#'
#' It is important to adhere to this prespecified order to enable \code{neImpute} to create valid pointers to these different types of predictor variables.
#' This requirement extends to the use of operators different than the \code{+} operator, such as the \code{:} and \code{*} operators (when e.g. adding interaction terms). 
#' For instance, the formula specifications \code{Y ~ X * M + C1 + C2}, \code{Y ~ X + M + X:M + C1 + C2} and \code{Y ~ X + X:M + M + C1 + C2} will create identical pointers to the different types of variables,
#' as the order of the unique predictor variables is identical in all three specifications. 
#' 
#' Furthermore, categorical exposures that are not coded as factors in the original dataset, should be specified as factors in the formula, 
#' using the \code{\link[base]{factor}} function, e.g. \code{Y ~ factor(X) + M + C1 + C2}. 
#' Quadratic or higher-order polynomial terms can be included as well, by making use of the \code{\link[base]{I}} function or by using the \code{\link[stats]{poly}} function.
#' For instance, \code{Y ~ X + I(X^2) + M + C1 + C2} and \code{Y ~ poly(X, 2, raw = TRUE) + M + C1 + C2} are equivalent and result in identical pointers to the different types of variables.
# We do not recommend the use of orthogonal polynomials (i.e. using the default argument specification \code{raw = FALSE} in \code{poly}).
#' 
#' The command \code{terms(object, "vartype")} (with \code{object} replaced by the name of the resulting expanded dataset) can be used to check whether valid pointers have been created.
#'
#' If multiple mediators are specified (\code{nMed > 1}), the natural indirect effect parameter in the natural effect model captures the joint mediated effect. That is, the effect of the exposure on the outcome via these mediators considered jointly. 
#' The remaining effect of the exposure on the outcome (not mediated through the specified mediators) is then captured by the natural indirect effect parameter.
#'
#' The type of imputation model can be defined by specifying an appropriate model-fitting function via the \code{FUN} argument (its default is \code{\link[stats]{glm}}).
#' This method can only be used with model-fitting functions that require a \code{formula} argument (so not when using e.g. \code{\link[SuperLearner]{SuperLearner}}).
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
#' @seealso \code{\link{neImpute}}, \code{\link{neImpute.default}}, \code{\link{neModel}}, \code{\link{expData}}
#' @examples
#' data(UPBdata)
#' 
#' ## example using glm imputation model with binary exposure
#' impData <- neImpute(UPB ~ attbin + negaff + gender + educ + age, 
#'                     family = binomial, data = UPBdata)
#' head(impData)
#' 
#' ## example using glm imputation model with continuous exposure
#' impData <- neImpute(UPB ~ att + negaff + gender + educ + age, 
#'                     family = binomial, data = UPBdata, nRep = 2)
#' head(impData)
#' 
#' ## example using vglm (yielding identical results as with glm)
#' library(VGAM)
#' impData2 <- neImpute(UPB ~ att + negaff + gender + educ + age, 
#'                      family = binomialff, data = UPBdata, 
#'                      nRep = 2, FUN = vglm)
#' head(impData2)
#' @export
neImpute.formula <- function (object, family, data, FUN = glm, nMed = 1, nRep = 5, 
    xSampling = c("quantiles", "random"), xFit, percLim = c(0.05, 
        0.95), ...) 
{
    args <- as.list(match.call())[-1L]
    if (missing(data)) 
      data <- environment(object)
    formula <- object
    if (missing(family)) 
        family <- formals(FUN)$family
    joint <- args$joint <- TRUE
    if (substitute(FUN) == "SuperLearner") {
        class(formula) <- "Yformula"
        Yind <- attr(neTerms(formula, nMed, joint = joint), "vartype")[["Y"]]
        Xind <- unlist(attr(neTerms(formula, nMed, joint = joint), 
            "vartype"))[-1]
        Y <- substitute(data[, Yind])
        X <- substitute(data[, Xind])
        args$object <- do.call(FUN, eval(list(Y = Y, X = X, family = family, 
            ...)))
        call <- substitute(list(Y = Y, X = X, family = family, 
            ...))
    }
    else {
        args$object <- do.call(FUN, eval(list(formula = formula, 
            family = family, data = data, ...)))
        call <- substitute(list(formula = formula, family = family, 
            data = data, ...))
        args$formula <- args$family <- args$data <- NULL
    }
    call[[1]] <- substitute(FUN)
    if (isS4(args$object)) 
        args$object@call <- call
    else args$object$call <- call
    expData <- do.call("neImpute.default", args)
    return(expData)
}
