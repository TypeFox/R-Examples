#' Linear hypotheses for natural effect models
#'
#' @description \code{neLht} allows to calculate linear combinations of natural effect model parameter estimates.\cr \code{neEffdecomp} automatically extracts relevant causal parameter estimates from a natural effect model.
#' @param model a fitted natural effect model object.
#' @param xRef a vector including reference levels for the exposure, \emph{x*} and \emph{x}, at which natural effect components need to be evaluated (see details). 
#' @param covLev a vector including covariate levels at which natural effect components need to be evaluated (see details). 
#' @param ... additional arguments (passed to \code{\link[multcomp]{glht}}).
#' @return An object of class \code{c("neLht", "glht")} (see \code{\link[multcomp]{glht}}). 
#' If the bootstrap is used for obtaining standard errors when fitting the \code{\link{neModel}} object, the returned object additionally inherits from class \code{"neLhtBoot"}. 
#' \code{neEffdecomp} returns an object that additionally inherits from class \code{"neEffdecomp"}.
#' 
#' See \code{\link{neLht-methods}} for methods for \code{neLht} objects (and \code{\link[=coef.glht]{glht-methods}} for additional methods for \code{glht} objects).
#' @details \code{neLht} is a wrapper of \code{\link[multcomp]{glht}} and offers the same functionality (see `Details' section of  \code{\link[multcomp]{glht}} for details on argument specification). 
#' It returns objects that inherit from the class \code{"neLht"} in order to make output of their corresponding methods (see \code{\link{neLht-methods}}) more compatible for natural effect models
#' containing bootstrap variance-covariance matrices and standard errors.
#' 
#' \code{neEffdecomp} is a convenience function that automatically extracts causal parameter estimates from a natural effect model
#' and derives natural effect components.
#' That is, natural direct, natural indirect and total causal effect estimates are returned if no exposure-mediator interaction is modelled (i.e. two-way decomposition). 
#' If mediated interaction is allowed for in the natural effect model, there are two ways of decomposing the total effect into (natural) direct and indirect effects components: 
#' either as the sum of the pure direct and the total indirect effect or as the sum of the pure indirect and the total direct effect (i.e. three-way decomposition).
#' In total, five causal effect estimates are returned in this case.
#' 
#' For continuous exposures, default exposure levels at which natural effect components are evaluated are \emph{x*} = 0 and \emph{x} = 1.
#' For multicategorical exposures, default levels are the reference level of the factor that encodes the exposure variable and the level corresponding to its first dummy variable for \emph{x*} and \emph{x}, respectively.  
#' If one wishes to evaluate natural effect components at different reference levels (e.g. if the natural effect model includes mediated interaction, quadratic or higher-order polynomial terms for the exposure; see examples), 
#' these can be specified as a vector of the form \code{c(x*,x)} via the \code{xRef} argument.
#' 
#' If applicable, covariate levels at which natural effect components are evaluated can also be specified. This is particularly useful for natural effect models encoding effect modification by baseline covariates (e.g. moderated mediation).
#' By default, these levels are set to 0 for continuous covariates and to the reference level for categorical covariates coded as factors. 
#' Different covariate levels can be specified via the \code{covLev} argument, which requires a vector including valid levels for covariates that are specified in the natural effect model (or a subset of covariates that are specified as modifiers of either the natural direct or indirect effect or both).
#' Levels need to be preceded by the name of the corresponding covariate, e.g., \code{covLev = c(gender = "M", age = 30)}. Covariates for which the levels are left unspecified are set to their default levels (see examples). 
#' The \code{\link{print}} and \code{\link[=summary.neLht]{summary}} functions for \code{neEffdecomp} objects return the covariate levels at which natural effect components are evaluated. 
#' No specific levels are returned for covariates that are not specified as modifier since effect decomposition is independent of the level of these covariates (see examples). 
#' 
#' @name neLht
#' @note \code{neEffdecomp} is internally called by \code{\link{plot.neModel}} to create confidence interval plots for \code{neModel} objects.
#' @seealso \code{\link{plot.neLht}}, \code{\link{neLht-methods}}, \code{\link[multcomp]{glht}}, \code{\link[=coef.glht]{glht-methods}}, \code{\link{neModel}}, \code{\link{plot.neModel}}, \code{\link{summary}}
#' @examples
#' data(UPBdata)
#' 
#' impData <- neImpute(UPB ~ att * negaff + gender + educ + age, 
#'                     family = binomial, data = UPBdata)
#' neMod <- neModel(UPB ~ att0 * att1 + gender + educ + age, 
#'                  family = binomial, expData = impData, se = "robust")
#' 
#' lht <- neLht(neMod, linfct = c("att0 = 0", "att0 + att0:att1 = 0", 
#'                                "att1 = 0", "att1 + att0:att1 = 0", 
#'                                "att0 + att1 + att0:att1 = 0"))
#' summary(lht)
#' 
#' ## or obtain directly via neEffdecomp
#' eff <- neEffdecomp(neMod)
#' summary(eff)
#' 
#' ## changing reference levels for multicategorical exposures
#' UPBdata$attcat <- factor(cut(UPBdata$att, 3), labels = c("L", "M", "H"))
#' impData <- neImpute(UPB ~ attcat * negaff + gender + educ + age,
#'                     family = binomial, data = UPBdata)
#' neMod <- neModel(UPB ~ attcat0 * attcat1 + gender + educ + age,
#'                  family = binomial, expData = impData, se = "robust")
#'                  
#' neEffdecomp(neMod)
#' neEffdecomp(neMod, xRef = c("L", "H"))
#' neEffdecomp(neMod, xRef = c("M", "H"))
#' 
#' 
#' ## changing reference levels for continuous exposures
#' impData <- neImpute(UPB ~ (att + I(att^2)) * negaff + gender + educ + age,
#'                     family = binomial, data = UPBdata)
#' neMod <- neModel(UPB ~ (att0 + I(att0^2)) * (att1 + I(att1^2)) + gender + educ + age,
#'                  family = binomial, expData = impData, se = "robust")
#' neEffdecomp(neMod)
#' neEffdecomp(neMod, xRef = c(-1, 0))
#' 
#' ## changing covariate levels when allowing for modification 
#' ## of the indirect effect by baseline covariates
#' impData <- neImpute(UPB ~ (att + negaff + gender + educ + age)^2,
#'                     family = binomial, data = UPBdata)
#' neMod <- neModel(UPB ~ att0 * att1 + gender + educ + age + att1:gender + att1:age,
#'                  family = binomial, expData = impData, se = "robust")
#' neEffdecomp(neMod)
#' neEffdecomp(neMod, covLev = c(gender = "F", age = 0)) # default covariate levels
#' neEffdecomp(neMod, covLev = c(gender = "M", age = 40))
#' neEffdecomp(neMod, covLev = c(gender = "M", age = 40, educ = "L"))
#' neEffdecomp(neMod, covLev = c(gender = "M", age = 40, educ = "M"))
#' neEffdecomp(neMod, covLev = c(gender = "M", age = 40, educ = "H"))
#' # effect decomposition is independent of education level
#' neEffdecomp(neMod, covLev = c(gender = "M")) 
#' # age is set to its default level when left unspecified
#' 
#' @export
NULL

#' Methods for linear hypotheses in natural effect models
#'
#' @description Obtain confidence intervals and statistical tests for linear hypotheses in natural effect models.
#' @param object an object of class \code{neLht}.
#' @param type the type of bootstrap intervals required. The default \code{"norm"} returns normal approximation bootstrap confidence intervals. Currently, \code{"norm"}, \code{"basic"}, \code{"perc"} and \code{"bca"} are supported (see \code{\link[boot]{boot.ci}}).
#' @param calpha a function computing the critical value. The default \code{univariate_calpha()} returns unadjusted confidence intervals, whereas \code{adjusted_calpha()} returns adjusted confidence intervals.
#' @param test a function for computing p-values. The default \code{univariate()} does not apply a multiple testing correction. The function \code{adjusted()} allows to correct for multiple testing (see \code{\link[multcomp]{summary.glht}} and \code{\link[multcomp]{adjusted}}) and \code{Chisquare()} allows to test global linear hypotheses.
#' @param ... additional arguments.
#' @inheritParams stats::confint.default
#' @name neLht-methods
#' @details 
#' \code{confint} yields bootstrap confidence intervals  or confidence intervals based on the sandwich estimator (depending on the type of standard errors requested when fitting the \code{\link{neModel}} object). 
#' Bootstrap confidence intervals are internally called via the \code{\link[boot]{boot.ci}} function from the \pkg{boot} package.
#' Confidence intervals based on the sandwich estimator are internally called via the corresponding \code{\link[multcomp]{confint.glht}} function from the \pkg{multcomp} package.
#' The default confidence level specified in \code{level} (which corresponds to the \code{conf} argument in \code{\link[boot]{boot.ci}}) is 0.95
#' and the default type of bootstrap confidence interval, \code{"norm"}, is based on the normal approximation.
#' Bias-corrected and accelerated (\code{"bca"}) bootstrap confidence intervals require a sufficiently large number of bootstrap replicates (for more details see \code{\link[boot]{boot.ci}}).
#'
#' A summary table with large sample tests, similar to that for \code{\link[multcomp]{glht}}, can be obtained using \code{summary}.
#' 
#' In contrast to \code{\link[multcomp]{summary.glht}}, which by default returns \emph{p}-values that are adjusted for multiple testing,
#' the summary function returns unadjusted \emph{p}-values. Adjusted \emph{p}-values can also be obtained by specifying the \code{test} argument 
#' (see \code{\link[multcomp]{adjusted}} for more details).
#' 
#' Global Wald tests considering all linear hypotheses simultaneously (i.e. testing the global null hypothesis) 
#' can be requested by specifying \code{test = Chisqtest()}.
#'
#' See \code{\link[=coef.glht]{glht-methods}} for additional methods for \code{glht} objects.
#'
#' @note For the bootstrap, \emph{z}-values in the summary table are simply calculated by dividing the parameter estimate by its corresponding bootstrap standard error. 
#' Corresponding \emph{p}-values in the summary table are only indicative, since the null distribution for each statistic is assumed to be approximately standard normal.
#' Therefore, whenever possible, it is recommended to focus mainly on bootstrap confidence intervals for inference, rather than the provided \emph{p}-values.
#' @seealso \code{\link{neLht}}, \code{\link{plot.neLht}}, \code{\link[multcomp]{glht}}, \code{\link[=coef.glht]{glht-methods}}
#' @examples
#' data(UPBdata)
#' 
#' impData <- neImpute(UPB ~ att * negaff + gender + educ + age, 
#'                     family = binomial, data = UPBdata)
#' neMod <- neModel(UPB ~ att0 * att1 + gender + educ + age, 
#'                  family = binomial, expData = impData, se = "robust")
#' 
#' lht <- neLht(neMod, linfct = c("att0 = 0", "att0 + att0:att1 = 0", 
#'                                "att1 = 0", "att1 + att0:att1 = 0", 
#'                                "att0 + att1 + att0:att1 = 0"))
#' 
#' ## obtain confidence intervals
#' confint(lht)
#' confint(lht, parm = c("att0", "att0 + att0:att1"))
#' confint(lht, parm = 1:2, level = 0.90)
#' 
#' ## summary table
#' summary(lht)
#' 
#' ## summary table with omnibus Chisquare test
#' summary(lht, test = Chisqtest())
NULL

#' Confidence interval plots for linear hypotheses in natural effect models
#'
#' @description Confidence interval plots for linear hypotheses in natural effect models.
#' @param x an object of class \code{neLht}.
#' @param ci.type the type of bootstrap intervals required (see \code{type} argument in \code{\link[medflex]{neModel-methods}}).
#' @param transf transformation function to be applied internally on the (linear hypothesis) estimates and their confidence intervals (e.g. \code{exp} for logit or Poisson regression). The default is \code{identity} (i.e. no transformation).
#' @param ylabels character vector containing the labels for the (linear hypothesis) estimates to be plotted on the y-axis.
#' @param yticks.at numeric vector containing the y-coordinates (from 0 to 1) to draw the tick marks for the different estimates and their corresponding confidence intervals.
#' @param ... additional arguments.
#' @inheritParams stats::confint.default
#' @name plot.neLht
#' @details This function is an adapted version of \code{\link[multcomp]{plot.glht}} from the \pkg{multcomp} package and
#' yields confidence interval plots for each of the linear hypothesis parameters.
#' @seealso \code{\link{neModel}}, \code{\link{neLht}}, \code{\link{neEffdecomp}}
#' @examples
#' data(UPBdata)
#' 
#' impData <- neImpute(UPB ~ att * negaff + gender + educ + age, 
#'                     family = binomial, data = UPBdata)
#' neMod <- neModel(UPB ~ att0 * att1 + gender + educ + age, 
#'                  family = binomial, expData = impData, se = "robust")
#' 
#' lht <- neLht(neMod, linfct = c("att0 = 0", "att0 + att0:att1 = 0", 
#'                                "att1 = 0", "att1 + att0:att1 = 0", 
#'                                "att0 + att1 + att0:att1 = 0"))
#' 
#' ## all pairs return identical output
#' plot(confint(lht), transf = exp)
#' plot(lht, transf = exp)
#' 
#' plot(neEffdecomp(neMod), transf = exp)
#' plot(neMod, transf = exp)
#' 
#' \dontshow{
#'   plot(neEffdecomp(neMod), level = 0.8, transf = exp, ylabels = c("PDE", "TDE", "PIE", "TIE", "TE"), yticks.at = c(0, 0.1, 0.5, 0.6, 1))
#'   plot(neMod, level = 0.8, transf = exp, ylabels = c("PDE", "TDE", "PIE", "TIE", "TE"), yticks.at = c(0, 0.1, 0.5, 0.6, 1))
#'   
#'   lht <- neLht(neMod, linfct = c("att0 = 0"))
#'   summary(lht)
#'   lht <- neLht(neMod, linfct = c("att0 = 0", "att0 + att0:att1 = 2"))
#'   summary(lht)
#' }
#' @export
NULL

coef.neLhtCI <- function (object, ...) 
{
    attr(object, "coef")
}

#' @rdname neLht-methods
#' @export
confint.neLhtBoot <- function (object, parm, level = 0.95, type = "norm", ...) 
{
    if (missing(parm)) 
        parm <- names(coef(object))
    else if (is.numeric(parm)) 
        parm <- names(coef(object))[parm]
    extractBootci <- function(index, x, level, type) rev(rev(boot::boot.ci(x, 
        conf = level, type = type, t0 = object$linfct[index, 
            ] %*% object$model$bootRes$t0, t = as.vector(object$linfct[index, 
            ] %*% t(object$model$bootRes$t)))[[4]])[1:2])
    ci <- t(sapply(seq.int(nrow(object$linfct)), extractBootci, 
        object$model$bootRes, level = level, type = type))
    dimnames(ci) <- list(rownames(object$linfct), paste0(100 * 
        level, c("% LCL", "% UCL")))
    ci <- ci[parm, ]
    attributes(ci) <- c(attributes(ci), list(level = level, coef = coef(object)[parm], 
        R = attr(object, "R"), type = type))
    class(ci) <- c("neBootCI", "neLhtCI", class(ci))
    return(ci)
}

#' @rdname neLht-methods
#' @export
confint.neLht <- function (object, parm, level = 0.95, calpha = univariate_calpha(), ...) 
{
    class(object) <- c("glht", class(object))
    ci <- confint(object, level = level, calpha = calpha, ...)$confint[, c("lwr", "upr")]
    dimnames(ci)[[2]] <- paste0(100 * level, c("% LCL", "% UCL"))
    ci <- ci[parm, ]
    attributes(ci) <- c(attributes(ci), list(level = level, coef = coef(object)[parm], calpha = calpha))
    class(ci) <- c("neLhtCI", class(ci))
    return(ci)
}

#' @rdname neLht
#' @export
neEffdecomp <- function (model, xRef, covLev, ...) 
{
    UseMethod("neEffdecomp")
}

#' @rdname neLht
#' @export
neEffdecomp.neModel <- function (model, xRef, covLev, ...) 
{
    xFact <- is.factor(model$neModelFit$data[, attr(terms(model), "vartype")$Xexp[[1]]])
    if (xFact) xRefCheck <- levels(model$neModelFit$data[, attr(terms(model), "vartype")$Xexp[[1]]])
    if (missing(xRef)) {
      xRef <- if (xFact) xRefCheck[1:2] else c(0, 1)
    } else {
      if (xFact && !all(xRef %in% xRefCheck)) {
        warning(gettextf("Invalid reference levels! Default reference levels %s were used instead.", 
                         paste0("c(", paste(paste0("'", xRefCheck, "'"), collapse = ", "), ")")))
        xRef <- xRefCheck[1:2]
      }
    }
    calcContr <- function(x, formula, covDat) {
      if (xFact) x <- factor(x, levels = xRefCheck)
      dat1 <- if (nrow(covDat)) data.frame(1, x[1], x[2], covDat) else data.frame(1, x[1], x[2])
      names(dat1) <- all.vars(formula)
      oldVars <- grep("factor\\(", dimnames(attr(terms(formula), "factors"))[[1]], value = TRUE)
      if (length(oldVars)) {
        newVars <- names(which(sapply(all.vars(formula), grep, oldVars)==TRUE))
        tmp <- mgsub(oldVars, newVars, deparse(formula[[3]]), fixed = TRUE)
        formula[[3]] <- as.call(parse(text = tmp))[[1]]
      } #
      modmat1 <- model.matrix(formula, data = dat1)
      dat2 <- if (nrow(covDat)) data.frame(1, x[3], x[4], covDat) else data.frame(1, x[3], x[4])
      names(dat2) <- all.vars(formula)
      modmat2 <- model.matrix(formula, data = dat2)
      return(t(modmat1 - modmat2))
    }
    list <- list(xRef[c(2, 1, 1, 1)],
                 xRef[c(2, 2, 1, 2)],
                 xRef[c(1, 2, 1, 1)],
                 xRef[c(2, 2, 2, 1)],
                 xRef[c(2, 2, 1, 1)])
    form <- model$neModelFit$formula
    vartype <- attr(terms(model), "vartype")
    if (is.na(vartype$Y)) vartype$Y <- as.character(form[[2]])
    cov <- all.vars(form)[!all.vars(form) %in% unlist(vartype[c("Y", "Xexp")])]
    varterms <- dimnames(attr(terms(form), "factors"))[[1]]
    covTerms <- varterms[sapply(cov, grep, varterms)]
    if (identical(cov, covTerms)) {
        dat <- substitute(model$neModelFit$data)
    } 
    else {
        modframe <- model.frame(model$neModelFit, data = model$neModelFit$data)
        dat <- substitute(modframe)
    }
    if (!missing(covLev) & !length(cov)) warning("As the natural effect model does not encode covariate effects, specified covariate levels are not taken into account.")
    tmp <- attr(terms(form), "factors")[covTerms, which(attr(terms(form), "order") > 1), drop = FALSE]
    covModifier <- if (!is.null(tmp)) cov[unlist(sapply(cov, grep, names(which(rowSums(as.matrix(tmp)) > 0))))] else NULL
    covClass <- sapply(eval(dat)[, covTerms, drop = FALSE], class)
    if (missing(covLev)) {
      covLev <- rep(NA, length(cov))
      names(covLev) <- cov
    }
    if (is.null(names(covLev))) {
      warning("Please provide names for the covariates! The specified covariate levels were discarded and default levels were used instead.")
    } 
    else {
      if (!all(names(covLev) %in% cov)) warning("For some covariate levels, the corresponding covariate name was either not provided or invalid. These levels were discarded and default levels were used instead.")
    }
    covLev <- covLev[cov]
    covMat <- cbind(covLev, cov, covClass, covTerms)
    covList <- split(covMat, 1:nrow(covMat))
    lapply(covList, function(x) {
      switch(x[3],
             "factor" = {if (all(!is.na(x[1]), !x[1] %in% levels(eval(dat)[, x[4]])))  
               warning(gettextf("Invalid covariate levels for %s! Default levels for this covariate were used instead.", x[2]))},
             "numeric" = {if (all(!is.na(x[1]), is.na(suppressWarnings(as.numeric(x[1])))))
               warning(gettextf("Invalid covariate levels for %s! Default levels for this covariate were used instead.", x[2]))},
             "integer" = {if (all(!is.na(x[1]), is.na(suppressWarnings(as.numeric(x[1])))))
               warning(gettextf("Invalid covariate levels for %s! Default levels for this covariate were used instead.", x[2]))})
    })
    covDat <- as.data.frame(lapply(covList, function(x) {
      switch(x[3],
             "factor" = {factor(if (all(!is.na(x[1]), x[1] %in% levels(eval(dat)[, x[4]]))) x[1] else levels(eval(dat)[, x[4]])[1], 
                                levels = levels(eval(dat)[, x[4]]))},
             "numeric" = {if (all(!is.na(x[1]), !is.na(as.numeric(x[1])))) as.numeric(x[1]) else 0},
             "integer" = {if (all(!is.na(x[1]), !is.na(as.numeric(x[1])))) as.numeric(x[1]) else 0})
    }))
    if (nrow(covDat)) {
      covLev <- as.matrix(covDat)
      colnames(covLev) <- cov
    }
    else {
      covLev <- NULL
    }
    K2 <- t(data.frame(lapply(list, calcContr, form, covDat)))
    K2 <- unique(K2)
    colnames(K2) <- names(coef(model)) #
    rownames <- if (nrow(K2) == 3) {
      c("natural direct effect", "natural indirect effect", "total effect")
    } else {
      c("pure direct effect", "total direct effect", "pure indirect effect", "total indirect effect", "total effect")
    }
    K <- matrix(0, nrow = nrow(K2), ncol = length(coef(model)), dimnames = list(rownames, names(coef(model))))
    K[, colnames(K2)] <- K2
    colnames(K) <- NULL
    effdecomp <- neLht(model, linfct = K)
    class(effdecomp) <- c("neEffdecomp", class(effdecomp))
    attributes(effdecomp) <- c(attributes(effdecomp), list(xRef = xRef, covLev = covLev, covModifier = covModifier))
    return(effdecomp)
}

#' @rdname neLht
#' @export
neLht <- function (model, ...) 
{
    UseMethod("neLht")
}

#' @rdname neLht
#' @export
neLht.neModel <- function (model, ...) 
{
    lht <- multcomp::glht(model, ...)
    if (inherits(model, "neModelBoot")) attr(lht, "R") <- model$bootRes$R
    class(lht) <- c(if (inherits(model, "neModelBoot")) "neLhtBoot", "neLht", class(lht))
    return(lht)
}

#' @rdname plot.neLht
#' @export
plot.neEffdecomp <- function (x, level = 0.95, transf = identity, 
    ylabels, yticks.at, ...) 
{
    args <- as.list(match.call())
    if (nrow(x$linfct) == 5)
        if (missing(yticks.at)) 
            args$yticks.at <- c(0, 0.15, 0.5, 0.65, 1)
        else if (missing(yticks.at)) 
            args$yticks.at <- NULL
    args[[1]] <- if (inherits(x, "neLhtBoot")) substitute(plot.neLhtBoot) else substitute(plot.neLht)
    eval(as.call(args))
}

#' @rdname plot.neLht
#' @export
plot.neLht <- function (x, level = 0.95, transf = identity, ylabels, 
    yticks.at, ...) 
{
    args <- as.list(match.call())[-1L]
    args <- c(substitute(confint), object = substitute(x),
              args[names(args) %in% names(formals(confint.neLhtBoot))])
    confint <- eval(as.call(args))
    plot(confint, transf = transf, ylabels = ylabels, yticks.at = yticks.at, 
        ...)
}

#' @rdname plot.neLht
#' @export
plot.neLhtBoot <- function (x, level = 0.95, ci.type = "norm", transf = identity, 
    ylabels, yticks.at, ...) 
{
    args <- as.list(match.call())[-1L]
    args <- c(substitute(confint), object = substitute(x), type = ci.type, 
              args[names(args) %in% names(formals(confint.neLhtBoot))])
    confint <- eval(as.call(args))
    plot(confint, transf = transf, ylabels = ylabels, yticks.at = yticks.at, 
         ...)
}

#' @export
plot.neLhtCI <- function (x, transf = identity, ylabels, yticks.at, main, xlab = NULL, 
    ylab = NULL, xlim, ylim, mar, mai, ...) 
{
    args <- as.list(match.call())[-1L]
    if (!grep("confint", args$x) == 1) {
        transf <- eval(args$x[[1]])
        args$x <- args$x[[-1]]
        x <- eval(args$x)
    }
    ci <- transf(x)
    xrange <- c(min(ci[, 1]), max(ci[, 2]))
    yvals <- yticks.at <- if (missing(yticks.at)) 
        nrow(ci):1
    else (nrow(ci) - 1) * (1 - yticks.at) + 1
    if (missing(ylabels)) 
        ylabels <- dimnames(ci)[[1]]
    old.par <- par(no.readonly = TRUE)
    if (all(missing(mai), missing(mar))) {
        mar <- old.par$mar
        mar[2] <- (max(strwidth(ylabels, "inch")) + 0.4) * (old.par$mar/old.par$mai)[1]
    }
    else if (missing(mar)) {
        mar <- mai * (old.par$mar/old.par$mai)[1]
    }
    par(mar = mar)
    null <- transf(0)
    tmp.range <- c(min(null, xrange[1]), max(null, xrange[2]))
    tmp.range <- mean(tmp.range) + c(-1, 1) * 0.6 * diff(tmp.range)
    if (missing(xlab)) 
        xlab <- ""
    if (missing(ylab)) 
        ylab <- ""
    if (missing(xlim)) 
        xlim <- tmp.range
    if (missing(ylim)) 
        ylim <- c(0.5, nrow(ci) + 0.5)
    if (missing(main)) 
        main <- paste0(100 * attr(x, "level"), if (inherits(x, "neBootCI")) "% bootstrap CIs" else "% sandwich CIs")
    plot(c(ci[, 1], ci[, 2]), rep.int(yvals, 2), type = "n", 
        main = main, axes = FALSE, xlab = xlab, ylab = ylab, 
        xlim = xlim, ylim = ylim, ...)
    axis(1, las = 1, ...)
    axis(2, at = yvals, labels = ylabels, las = 1, ...)
    abline(v = null, lty = 2, lwd = 1, ...)
    left <- ci[, 1]
    left[!is.finite(left)] <- min(c(0, xlim[1] * 2))
    right <- ci[, 2]
    right[!is.finite(right)] <- max(c(0, xlim[2] * 2))
    segments(left, yvals, right, yvals, ...)
    points(transf(coef(x)), yvals, pch = 20, ...)
    box()
    par(mar = old.par$mar, mai = old.par$mai)
}

#' @method print neLht
#' @export
print.neLht <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    if (inherits(x, "neEffdecomp")) {
      cat("Effect decomposition on the scale of the linear predictor\n---\n")
      covModifier <- dimnames(attr(x, "covLev"))[[2]] %in% attr(x, "covModifier")
      sep <- if (all(covModifier) | !any(covModifier)) "" else ", "
      if (!is.null(attr(x, "covLev"))) cat("conditional on:", paste(paste(dimnames(attr(x, "covLev")[, covModifier, drop = FALSE])[[2]], attr(x, "covLev")[, covModifier, drop = FALSE], sep = " = ", collapse = ", "), 
                                                                    paste(dimnames(attr(x, "covLev")[, !covModifier, drop = FALSE])[[2]], collapse = ", "), sep = sep), "\n") 
      cat(paste0("with x* = ", attr(x, "xRef")[1], ", x = ", attr(x, "xRef")[2]), "\n---\n")
    }
    else cat("Linear hypotheses for natural effect models\n---\n")
    coef <- as.matrix(coef(x))
    dimnames(coef)[[2]] <- "Estimate"
    print(coef, digits = digits, ...)
}

#' @method print neLhtCI
#' @export
print.neLhtCI <- function (x, ...) 
{
  adj <- attr(attr(x, "calpha"), "type")
  attributes(x)[c("level", "coef", "class", "calpha")] <- NULL
  print.default(x)
  switch(adj, 
         "univariate" = cat("---\nconfidence intervals\nbased on the sandwich estimator\n\n"),
         "adjusted" = cat("---\nadjusted confidence intervals\nbased on the sandwich estimator\n\n"))
}

#' @method print summary.neLht
#' @export
print.summary.neLht <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    catSE <- if (inherits(x, "summary.neLhtBoot")) "with standard errors based on the non-parametric bootstrap\n---\n"
    else "with standard errors based on the sandwich estimator\n---\n"
    if (inherits(x$test, "mtest")) {
        if (inherits(x, "summary.neEffdecomp")) {
          cat("Effect decomposition on the scale of the linear predictor\n", catSE, sep = "")
          covModifier <- dimnames(attr(x, "covLev"))[[2]] %in% attr(x, "covModifier")
          sep <- if (all(covModifier) | !any(covModifier)) "" else ", "
          if (!is.null(attr(x, "covLev"))) cat("conditional on:", paste(paste(dimnames(attr(x, "covLev")[, covModifier, drop = FALSE])[[2]], attr(x, "covLev")[, covModifier, drop = FALSE], sep = " = ", collapse = ", "), 
                                                                        paste(dimnames(attr(x, "covLev")[, !covModifier, drop = FALSE])[[2]], collapse = ", "), sep = sep), "\n") 
          cat(paste0("with x* = ", attr(x, "xRef")[1], ", x = ", attr(x, "xRef")[2]), "\n---\n")  
        }
        else {
          cat("Linear hypotheses for natural effect models\n", catSE, sep = "")
        }
        if (!identical(x$rhs, rep(0, length(x$rhs)))) 
            dimnames(x$coef.table)[[1]] <- paste(dimnames(x$coef.table)[[1]], 
                "=", x$rhs)
        printCoefmat(x$coef.table, digits = digits, has.Pvalue = TRUE, 
            P.values = TRUE)
        switch(x$test$type, univariate = cat("(Univariate p-values reported)"), 
            `single-step` = cat("(Adjusted p-values reported -- single-step method)"), 
            Shaffer = cat("(Adjusted p-values reported -- Shaffer method)"), 
            Westfall = cat("(Adjusted p-values reported -- Westfall method)"), 
            cat("(Adjusted p-values reported --", x$test$type, 
                "method)"))
        cat("\n\n")
    }
    else if (inherits(x$test, "gtest")) {
        cat("Global linear hypothesis test for natural effect models\n", 
            catSE, sep = "")
        pr <- data.frame(x$test$SSH, x$test$df[1], x$test$pval)
        names(pr) <- c("Chisq", "DF", "Pr(>Chisq)")
        print(pr, digits = digits, ...)
    }
}

#' @rdname neLht-methods
#' @method summary neLht
#' @export
summary.neLht <- function (object, test = univariate(), ...) 
{
    class(object) <- c("glht", class(object))
    summary <- summary(object, test = test, ...)
    pq <- summary$test
    if (inherits(pq, "mtest")) {
        coef.table <- cbind(pq$coefficients, pq$sigma, pq$tstat, 
            pq$pvalues)
        dimnames(coef.table) <- list(names(coef(object)), c("Estimate", 
            "Std. Error", "z value", "Pr(>|z|)"))
        summary$coef.table <- coef.table
    }
    class(summary) <- c(if (inherits(object, "neEffdecomp")) "summary.neEffdecomp",
                        if (inherits(object, "neLhtBoot")) "summary.neLhtBoot", "summary.neLht", class(summary))
    return(summary)
}
