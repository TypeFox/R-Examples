#' Get summary information
#'
#' List the estimates, standard errors, p-values and confidence intervals for the coefficients of
#' robust linear regression models with compositional data as returned by \code{\link{complmrob}} or
#' \code{\link{bootcoefs}}
#'
#' @param object the object for which the summary information should be returned.
#' @param conf.level the level of the returned confidence intervals.
#' @param conf.type the type of the returned confidence interval (see \code{\link[boot]{boot.ci}} for the
#' meaning of this parameter).
#' @param ... ignored.
#'
#' @name summary-methods
NULL

# @describeIn summary for robust linear regression models with compositional data
#' @import robustbase
#' @rdname summary-methods
#' @export
summary.complmrob <- function(object, conf.level = 0.95, ...) {
    ret <- list(
        obj = object,
        ci = NULL,
        stats = NULL,
        r.squared = NULL,
        adj.r.squared = NULL,
        scale = NULL,
        type = "theoretical"
    );

    intercSe <- sqrt(object$models[[1]]$cov[object$coefind, object$coefind]);
    intercTval <- object$models[[1]]$coefficients[1] / intercSe;

    thParams <- list();

    if(object$intercept == TRUE) {
        thParams <- list("(Intercept)" = c(
            se = intercSe,
            tval = intercTval,
            pval = 2 * pt(abs(intercTval), object$models[[1]]$df.residual, lower.tail = FALSE)
        ));
    }

    thParams <- c(thParams, lapply(object$models, function(m) {
        se <- sqrt(m$cov[object$coefind, object$coefind]);
        tval <- m$coefficients[object$coefind] / se;
        return(c(
            se = se,
            tval <- tval,
            pval = 2 * pt(abs(tval), m$df.residual, lower.tail = FALSE)
        ));
    }));

    ret$stats <- as.data.frame(do.call(rbind, thParams))
    colnames(ret$stats) <- c("Std. Error", "t value", "Pr(>|t|)");

    ret$ci <- list();

    if(object$intercept == TRUE) {
        ret$ci <- list("(Intercept)" = confint(object$models[[1]], level = conf.level)[1L, ]);
    }

    ret$ci <- c(ret$ci, lapply(object$models, function(m) {
        confint(m, level = conf.level)[object$coefind, ]
    }));

    ret$ci <- do.call(rbind, ret$ci);

    sm <- summary(object$models[[1]]);
    ret$r.squared <- sm$r.squared;
    ret$adj.r.squared <- sm$adj.r.squared;
    ret$scale = object$models[[1]]$scale;

    class(ret) <- "summary.complmrob";
    return(ret);
}

# @describeIn summary for bootstrapped robust linear regression models with compositional data
#' @import robustbase
#' @rdname summary-methods
#' @export
summary.bccomplmrob <- function(object, conf.level = 0.95, conf.type = "perc", ...) {
    ret <- list(
        obj = object$model,
        R = object$R,
        ci = NULL,
        stats = NULL,
        r.squared = NULL,
        adj.r.squared = NULL,
        scale = NULL,
        type = "bootstrapped"
    );

    statBs <- do.call(rbind, lapply(object$bootres, function(bo) {
        bias <- mean(bo$t, na.rm = TRUE) - bo$t0;
        pval <- (sum((sign(bo$t0) * bo$t) < 0, na.rm = TRUE) + 1) / (sum(!is.na(bo$t)) + 1);
        c(bias = bias, se = sd(bo$t, na.rm = TRUE), pval = pval)
    }));

    ret$stats <- cbind("bias" = statBs[ , 1L], "Std. Error" = statBs[ , 2L], "Pr(b<>0)" = statBs[ , 3L]);

    ret$ci <- confint(object, level = conf.level, type = conf.type);

    sm <- summary(object$model$models[[1]]);
    ret$r.squared <- sm$r.squared;
    ret$adj.r.squared <- sm$adj.r.squared;
    ret$scale = object$model$models[[1]]$scale;

    class(ret) <- "summary.complmrob";
    return(ret);
}

# @describeIn summary for bootstrapped robust linear regression models
#' @import robustbase
#' @importFrom boot boot.ci
#' @rdname summary-methods
#' @export
summary.bclmrob <- function(object, conf.level = 0.95, conf.type = "perc", ...) {
    ret <- list(
        obj = object$model,
        R = object$R,
        ci = NULL,
        stats = NULL,
        r.squared = NULL,
        adj.r.squared = NULL,
        scale = NULL,
        type = "bootlmrob"
    );

    nc <- ncol(object$bootres$t);
    bsl <- split(object$bootres$t, rep.int(seq_len(nc), times = rep.int(nrow(object$bootres$t), nc)));

    statBs <- do.call(rbind, mapply(function(t, t0) {
        bias <- mean(t, na.rm = TRUE) - t0;
        pval <- (sum((sign(t0) * t) < 0, na.rm = TRUE) + 1) / (sum(!is.na(t)) + 1);
        c(bias = bias, se = sd(t, na.rm = TRUE), pval = pval)
    }, bsl, object$bootres$t0, SIMPLIFY = FALSE));

    ret$stats <- cbind("bias" = statBs[ , 1L], "Std. Error" = statBs[ , 2L], "Pr(b<>0)" = statBs[ , 3L]);

    ret$ci <- confint(object, level = conf.level, type = conf.type);

    sm <- summary(object$model);
    ret$r.squared <- sm$r.squared;
    ret$adj.r.squared <- sm$adj.r.squared;
    ret$scale = object$model$scale;

    class(ret) <- "summary.complmrob";
    return(ret);
}

#' Print the summary information
#'
#' @param x the summary object.
#' @param digits the number of digits for the reported figures
#' @param signif.stars should stars be displayed to show the significance of certain figures
#' @param ... further arguments currently not used
#' @export
print.summary.complmrob <- function(x, digits = max(3, getOption("digits") - 3), signif.stars = getOption("show.signif.stars"), ...) {
    cat("Robust linear regression with compositional covariates\n");
    print(x$obj$call);
    if(x$type == "bootstrapped" || x$type == "bootlmrob") {
        cat("\nStandard errors and derived statistics are base on ", x$R, " bootstrap replications\n", sep = "");
    }
    cat("\nCoefficients:\n");

    prdf <- cbind("Estimate" = coef(x$obj), x$stats);
    printCoefmat(prdf, digits = digits, signif.stars = signif.stars, cs.ind = seq_len(ncol(prdf) - 1), tst.ind = NULL)

    cat("\nConfidence intervals:\n");
    printCoefmat(cbind("Estimate" = coef(x$obj), x$ci), P.values = FALSE, cs.ind = seq_len(3), tst.ind = NULL);

    cat("\nRobust residual standard error:", format(signif(x$scale, digits)), "\n");
    cat("Multiple R-squared: ", formatC(x$r.squared, digits = digits), sep = "");
    cat("\tAdjusted R-squared: ", formatC(x$adj.r.squared, digits = digits), "\n", sep = "");
}
