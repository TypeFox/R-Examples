################################################################################
##  Author: Benjamin Hofner, benjamin.hofner@fau.de

################################################################################
# Prettify function for summary tables
prettify <- function(object, ...)
    UseMethod("prettify")

prettify.summary.lm <- function(object, labels = NULL, sep = ": ", extra.column = FALSE,
                                confint = TRUE, level = 0.95,
                                smallest.pval = 0.001, digits = NULL, scientific = FALSE,
                                signif.stars = getOption("show.signif.stars"), ...) {

    .call <- match.call()
    res <- as.data.frame(coef(object))

    ## compute confidence interval or extract it from confint
    if (is.logical(confint)) {
        if (confint) {
            mod <- refit_model(cl = object$call,
                               ENV = attr(object$terms, ".Environment"),
                               summary = object, .call = .call)
            if (is.logical(mod)) {
                ## model could not be refitted, i.e., mod == FALSE
                confint <- mod
            } else {
                CI <- confint(mod, level = level)
            }
        }
    } else {
        CI <- confint
        confint <- TRUE
    }

    if (confint){
        res$CI_lower <- CI[,1]
        res$CI_upper <- CI[,2]
        ## move confint to the front
        newVars <- (ncol(res) -1):ncol(res)
        res <- cbind(res[, 1, drop = FALSE],
                     res[, newVars],
                     res[, - c(1, newVars)])
        names(res)[2] <- "CI (lower)"
        names(res)[3] <- "CI (upper)"
    }

    ## use variable names as labels
    if (is.null(labels)) {
        labels <- names(attr(object$terms, "dataClasses"))
        names(labels) <- labels
    }

    prettify(res, labels = labels, sep = sep, extra.column = extra.column,
             smallest.pval = smallest.pval, digits = digits,
             scientific = scientific, signif.stars = signif.stars, ...)
}

prettify.summary.glm <- function(object, labels = NULL, sep = ": ", extra.column = FALSE,
                                 confint = TRUE, level = 0.95, OR = TRUE,
                                 smallest.pval = 0.001, digits = NULL, scientific = FALSE,
                                 signif.stars = getOption("show.signif.stars"),
                                 ...) {

    .call <- match.call()
    res <- as.data.frame(coef(object))
    if (OR <- (object$family$family == "binomial" && OR)) {
        res$OR <- exp(res$Estimate)
    }

    ## compute confidence interval or extract it from confint
    if (is.logical(confint)) {
        if (confint) {
            mod <- refit_model(cl = object$call,
                               ENV = attr(object$terms, ".Environment"),
                               summary = object, .call = .call)
            if (is.logical(mod)) {
                ## model could not be refitted, i.e., mod == FALSE
                confint <- mod
            } else {
                CI <- confint(mod, level = level)
            }
        }
    } else {
        CI <- confint
        confint <- TRUE
    }

    if (confint){
        if (OR) {
            res$CI_lower <- exp(CI[,1])
            res$CI_upper <- exp(CI[,2])
            ## move confint to the front
            newVars <- (ncol(res) - 2):ncol(res)
            res <- cbind(res[, 1, drop = FALSE],
                         res[, newVars],
                         res[, - c(1, newVars)])
            names(res)[2] <- "Odds Ratio"
            names(res)[3] <- "CI (lower)"
            names(res)[4] <- "CI (upper)"
        } else {
            res$CI_lower <- CI[,1]
            res$CI_upper <- CI[,2]
            ## move confint to the front
            newVars <- (ncol(res) -1):ncol(res)
            res <- cbind(res[, 1, drop = FALSE],
                         res[, newVars],
                         res[, - c(1, newVars)])
            names(res)[2] <- "CI (lower)"
            names(res)[3] <- "CI (upper)"
        }
    }

    ## use variable names as labels
    if (is.null(labels)) {
        labels <- names(attr(object$terms, "dataClasses"))
        names(labels) <- labels
    }

    prettify(res, labels = labels, sep = sep, extra.column = extra.column,
             smallest.pval = smallest.pval, digits = digits,
             scientific = scientific, signif.stars = signif.stars, ...)
}

prettify.summary.coxph.penal <- prettify.summary.coxph <-
    function(object, labels = NULL, sep = ": ", extra.column = FALSE,
             confint = TRUE, level = 0.95, HR = TRUE,
             smallest.pval = 0.001, digits = NULL, scientific = FALSE,
             signif.stars = getOption("show.signif.stars"),
             env = parent.frame(), ...) {

    .call <- match.call()
    res <- as.data.frame(coef(object))
    if (!HR) {
        res$"exp(coef)" <- NULL
    } else {
        if (is.null(res$"exp(coef)"))
            res$"exp(coef)" <- exp(res$coef)
    }
    if (is.null(labels) || (is.logical(confint) && confint)) {
        mod <- refit_model(cl = object$call, ENV = env,
                           summary = object, .call = .call)
    }
    if (is.null(labels) && is.logical(mod))
        stop("Model can't be refitted and no labels are specified. ",
             "Please specify labels.")


    ## compute confidence interval or extract it from confint
    if (is.logical(confint)) {
        if (confint) {
            if (is.logical(mod)) {
                ## model could not be refitted, i.e., mod == FALSE
                confint <- mod
            } else {
                CI <- confint(mod, level = level)
            }
        }
    } else {
        CI <- confint
        confint <- TRUE
    }

    if (confint){
        message("Confidence intervals are experimental only;\n",
                "Model refitted but original environment not available.\n")
        res$CI_upper <- res$CI_lower <- NA
        if (HR) {
            res$CI_lower[1:nrow(CI)] <- exp(CI[,1])
            res$CI_upper[1:nrow(CI)] <- exp(CI[,2])
            ## move confint to the front
            res <- cbind(res[, c("coef", "exp(coef)"), drop = FALSE],
                         res[, c("CI_lower", "CI_upper")],
                         res[, !colnames(res) %in% c("coef", "exp(coef)", "CI_lower", "CI_upper", "se2")])
            names(res)[2] <- "Hazard Ratio"
            names(res)[3] <- "CI (lower)"
            names(res)[4] <- "CI (upper)"
        } else {
            res$CI_lower[1:nrow(CI)] <- CI[,1]
            res$CI_upper[1:nrow(CI)] <- CI[,2]
            ## move confint to the front
            res <- cbind(res[, c("coef"), drop = FALSE],
                         res[, c("CI_lower", "CI_upper")],
                         res[, !colnames(res) %in% c("coef", "CI_lower", "CI_upper", "se2")])
            names(res)[2] <- "CI (lower)"
            names(res)[3] <- "CI (upper)"
        }
    }

    ## use variable names as labels
    if (is.null(labels)) {
        labels <- names(attr(mod$terms, "dataClasses"))
        names(labels) <- labels
    }

    prettify(res, labels = labels, sep = sep, extra.column = extra.column,
             smallest.pval = smallest.pval, digits = digits,
             scientific = scientific, signif.stars = signif.stars, ...)
}

prettify.summary.lme <- function(object, labels = NULL, sep = ": ", extra.column = FALSE,
                                 confint = TRUE, level = 0.95,
                                 smallest.pval = 0.001, digits = NULL, scientific = FALSE,
                                 signif.stars = getOption("show.signif.stars"),
                                 ...) {

    .call <- match.call()
    res <- as.data.frame(object$tTable)

    ## compute confidence interval or extract it from confint
    if (is.logical(confint)) {
        if (confint) {
            mod <- refit_model(cl = object$call,
                               ENV = attr(object$terms, ".Environment"),
                               summary = object, .call = .call)
            if (is.logical(mod)) {
                ## model could not be refitted, i.e., mod == FALSE
                confint <- mod
            } else {
                CI <- confint(mod, level = level)
            }
        }
    } else {
        CI <- confint
        confint <- TRUE
    }

    if (confint){
        res$CI_lower <- CI[,1]
        res$CI_upper <- CI[,2]
        ## move confint to the front
        newVars <- (ncol(res) -1):ncol(res)
        res <- cbind(res[, 1, drop = FALSE],
                     res[, newVars],
                     res[, - c(1, newVars)])
        names(res)[2] <- "CI (lower)"
        names(res)[3] <- "CI (upper)"
    }

    ## use variable names as labels
    if (is.null(labels)) {
        labels <- names(attr(object$terms, "dataClasses"))
        names(labels) <- labels
    }

    prettify(res, labels = labels, sep = sep, extra.column = extra.column,
             smallest.pval = smallest.pval, digits = digits,
             scientific = scientific, signif.stars = signif.stars, ...)
}

prettify.summary.merMod <- function(object,
                     labels = NULL, sep = ": ", extra.column = FALSE,
                     confint = TRUE, level = 0.95,
                     smallest.pval = 0.001, digits = NULL, scientific = FALSE,
                     signif.stars = getOption("show.signif.stars"),
                     method = c("profile", "Wald", "boot"), B = 1000,
                     env = parent.frame(), ...) {

    .call <- match.call()
    res <- as.data.frame(coefficients(object))

    if (is.null(labels) || (is.logical(confint) && confint)) {
        mod <- refit_model(cl = object$call, ENV = env,
                           summary = object, .call = .call)
    }
    if (is.null(labels) && is.logical(mod))
        stop("Model can't be refitted and no labels are specified. ",
             "Please specify labels.")

    ## compute confidence interval or extract it from confint
    if (is.logical(confint)) {
        if (confint) {
            if (is.logical(mod)) {
                ## model could not be refitted, i.e., mod == FALSE
                confint <- mod
            } else {
                CI <- confint(mod, level = level, method = method, nsim = B,
                              ...)[rownames(res), ]
            }
        }
    } else {
        CI <- confint
        confint <- TRUE
    }

    if (confint){
        message("Confidence intervals are experimental only;\n",
                "Model refitted but original environment not available.\n")
        res$CI_lower <- CI[,1]
        res$CI_upper <- CI[,2]
        ## move confint to the front
        newVars <- (ncol(res) -1):ncol(res)
        res <- cbind(res[, 1, drop = FALSE],
                     res[, newVars],
                     res[, - c(1, newVars)])
        names(res)[2] <- "CI (lower)"
        names(res)[3] <- "CI (upper)"
    }

    ## use variable names as labels
    if (is.null(labels)) {
        labels <- names(attr(attr(mod@frame, "terms"), "dataClasses"))
        names(labels) <- labels
    }

    prettify(res, labels = labels, sep = sep, extra.column = extra.column,
             smallest.pval = smallest.pval, digits = digits,
             scientific = scientific, signif.stars = signif.stars, ...)
}

## nocov start (exclude this function from test coverage)
## function for lme4 version < 1.0 only
prettify.summary.mer <- function(object,
                     labels = NULL, sep = ": ", extra.column = FALSE,
                     confint = TRUE, level = 0.95,
                     smallest.pval = 0.001, digits = NULL, scientific = FALSE,
                     signif.stars = getOption("show.signif.stars"),
                     simulate = c("ifneeded", TRUE, FALSE), B = 1000, ...) {

    .call <- match.call()
    res <- as.data.frame(object@coefs)

    ## compute confidence interval or extract it from confint
    if (is.logical(confint)) {
        if (confint) {
            mod <- refit_model(cl = object@call,
                               ENV = attr(attr(object@frame, "terms"), ".Environment"),
                               summary = object, .call = .call)
            if (is.logical(mod)) {
                ## model could not be refitted, i.e., mod == FALSE
                confint <- mod
            } else {
                CI <- confint(mod, level = level, simulate = simulate, B = B, ...)
            }
        }
    } else {
        CI <- confint
        confint <- TRUE
    }

    if (confint){
        res$CI_lower <- CI[,1]
        res$CI_upper <- CI[,2]
        ## move confint to the front
        newVars <- (ncol(res) -1):ncol(res)
        res <- cbind(res[, 1, drop = FALSE],
                     res[, newVars],
                     res[, - c(1, newVars)])
        names(res)[2] <- "CI (lower)"
        names(res)[3] <- "CI (upper)"
    }

    ## use variable names as labels
    if (is.null(labels)) {
        labels <- names(attr(attr(object@frame, "terms"), "dataClasses"))
        names(labels) <- labels
    }

    prettify(res, labels = labels, sep = sep, extra.column = extra.column,
             smallest.pval = smallest.pval, digits = digits,
             scientific = scientific, signif.stars = signif.stars, ...)
}
## nocov end

prettify.anova <- function(object, labels = NULL,
                           smallest.pval = 0.001, digits = NULL, scientific = FALSE,
                           signif.stars = getOption("show.signif.stars"), ...){

    res <- as.data.frame(object)
    res <- prettifyPValue(res, smallest.pval, digits, scientific, signif.stars, ...)

    res$varlabel <- dimnames(object)[[1]]
    res$varlabel <- as.character(res$varlabel)
    newVars <- ncol(res)
    res <- cbind(res[, newVars],
                 res[, - newVars])
    names(res)[1] <- " "
    rownames(res) <- NULL

    if (!is.null(labels)) {
        idx <- res[, 1] %in% names(labels)
        if (any(idx == TRUE))
            res[, 1] <- as.character(res[, 1])
        res[idx, 1] <- labels[res[idx, 1]]
    }

    res <- res[!apply(res, 1, function(x) any(is.na(x))), ]
    res
}

prettify.data.frame <- function(object, labels = NULL, sep = ": ", extra.column = FALSE,
                                smallest.pval = 0.001, digits = NULL, scientific = FALSE,
                                signif.stars = getOption("show.signif.stars"),
                                ...) {
    ## get row names
    nms <- new_nms <- rownames(object)

    if (is.null(labels)) {
        if (extra.column)
            warning(sQuote("extra.column"),
                    " can only be used if labels are specified")
        extra.column <- FALSE
    } else {
        ## order labels to avoid matching with substrings
        labels <- labels[rev(order(sapply(names(labels), nchar)))]
    }

    ## make extra column for factor levels if needed
    if (extra.column) {
        object$varlabel <- " "
        object$"FactorLevel" <- " "
        ## move Factor Levels to the front
        newVars <- (ncol(object) -1):ncol(object)
        object <- cbind(object[, newVars],
                        object[, - newVars])
        names(object)[1] <- " "
        object[,1] <- as.character(object[,1])
        names(object)[2] <- "Factor Level"
        object[,2] <- as.character(object[,2])
    } else {
        object$varlabel <- new_nms
        newVars <- ncol(object)
        object <- cbind(object[, newVars],
                        object[, - newVars])
        names(object)[1] <- " "
        object[,1] <- as.character(object[,1])
    }

    if (!is.null(labels)) {
        for (i in 1:length(labels)) {
            idx <- grep(names(labels)[i], nms)
            if (!length(idx) == 0){
                ## Is there a factor level?
                if (any(grepl(paste("^",names(labels)[i], "$", sep = ""), nms[idx]))) {
                    ## if not replace variable names with labels
                    new_nms[idx] <- gsub(names(labels)[i], labels[i], nms[idx])
                } else {
                    ## if factors are present separate variable name and factor
                    ## level
                    if (extra.column) {
                        ## replace variable name with label and discard
                        ## everything else
                        new_nms[idx] <- gsub(paste("^",names(labels)[i], "(.*)", sep = ""),
                                             labels[i],
                                             nms[idx])
                        ## remove duplicate variable labels
                        new_nms[idx][duplicated(new_nms[idx])] <- ""
                        ## extract variable levels
                        object[idx, 2] <- gsub(paste("^",names(labels)[i], "(.*)", sep = ""),
                                               "\\1",
                                               nms[idx])
                    } else {
                        new_nms[idx] <- gsub(paste("^",names(labels)[i], "(.*)", sep = ""),
                                             paste(labels[i], sep, "\\1", sep = ""),
                                             nms[idx])
                    }
                }
                nms[idx] <- ""
            }
        }
    }
    object[, 1] <- new_nms
    rownames(object) <- NULL

    object <- prettifyPValue(object, smallest.pval, digits, scientific,
                             signif.stars, ...)

    object

}


### helper for pretty p-values
prettifyPValue <- function(object,
                           smallest.pval = 0.001, digits = NULL, scientific = FALSE,
                           signif.stars = getOption("show.signif.stars"), ...) {

    wchPval <- grep("Pr(.*)|p-value|^p$", names(object))
    if (length(wchPval) != 0) {
        if (signif.stars) {
            object$signif <- symnum(object[, wchPval], corr = FALSE, na = FALSE,
                                    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                                    symbols = c("***", "**", "*", ".", " "))
            names(object)[names(object) == "signif"] <- "   "
        }
        r.digits <- 10
        num <- strsplit(as.character(smallest.pval), "\\.")[[1]]
        if (!is.null(num[2]))
            r.digits <- nchar(num[2])
        object[, wchPval] <- format.pval(round(object[, wchPval], digits = r.digits),
                                         digits = digits, scientific = scientific,
                                         eps = smallest.pval, ...)
    }

    if (!is.null(digits)) {
        for (i in names(object)[-wchPval]) {
            if (is.numeric(object[, i]))
                object[, i] <- format(zapsmall(object[, i]), digits = digits,
                                      scientific = scientific, ...)
        }
    }

    object
}
