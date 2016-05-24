#' Extraction of p-value from a statistical test
#'
#' These calculations are written in such a way that they avoid rounding off errors
#' that plague the \pkg{survival} and \pkg{cmprsk} packages. 
#'
#' @param x Test, i.e. a fitted object of a supported type.
#' @param log_p Whether to return the logarithm of the p-value.
#' @param ... Sent to class method.
#' @return p-value.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @seealso \code{\link{pvalue.crr}}, \code{\link{pvalue.survdiff}},
#'   \code{\link{pvalue.cuminc}}
#' @export
pvalue <- function(x, log_p=FALSE, ...) UseMethod("pvalue")


#' Extract p-value from a Cox proportional hazards model
#' 
#' Based on \code{\link{summary.coxph}}.
#'
#' @method pvalue coxph
#' @param x Fitted \code{\link[survival]{coxph}} model.
#' @param log_p Whether to return the logarithm of the p-value.
#' @param test What test to calculate. \code{"likelihood"} is short for means
#'   likelihood ratio test.
#' @param ... Ignored. Kept for S3 consistency.
#' @return p-value.
#' @seealso \code{\link{pvalue}}
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
pvalue.coxph <- function(x, log_p=FALSE, test=c("logrank", "wald", "likelihood"), ...){
    test <- match.arg(test)
    df <- sum(!is.na(x$coefficients))
    switch(test,
        logrank = pchisq(x$score, df, lower.tail=FALSE, log.p=log_p),
        wald = pchisq(x$score, df, lower.tail=FALSE, log.p=log_p),
        likelihood = pchisq(x$logtest, df, lower.tail=FALSE, log.p=log_p))
}


#' Extract p-value from a cumulative incidence estimation
#' 
#' This is also known as Gray's test.
#' 
#' @method pvalue cuminc
#' @param x Fitted \code{\link[cmprsk]{cuminc}} estimate.
#' @param log_p Whether to return the logarithm of the p-value.
#' @param ... Ignored. Kept for S3 consistency.
#' @return p-value.
#' @seealso \code{\link{pvalue}}
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
pvalue.cuminc <- function(x, log_p=FALSE, ...){
    pchisq(x$Tests[,"stat"], x$Tests[,"df"], lower.tail=FALSE, log.p=log_p)
}


#' Extracts p-value from a competing risk model
#' 
#' @method pvalue crr
#' @param x Fitted crr model, as returned by \code{\link[cmprsk]{crr}}.
#' @param log_p Whether to return the logarithm of the p-value.
#' @param ... Ignored. Kept for S3 consistency.
#' @return Two-sided p-value.
#' @examples
#' library(cmprsk)
#' time <- 1:20
#' event <- c(rep(0, 9), rep(2, 3), rep(1, 8))
#' data <- rep(0:1, each=10)
#' x <- crr(time, event, data)
#' 
#' # Compare p-values of implementations
#' print(x)
#' pvalue(x)
#' @seealso \code{\link{pvalue}}
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
pvalue.crr <- function(x, log_p=FALSE, ...){
    pval <- pnorm(abs(x$coef)/sqrt(diag(x$var)), lower.tail=FALSE, log.p=log_p)
    if(log_p) log(2) + pval else 2 * pval
}

#' Extracts p-value from a logrank test
#' 
#' @method pvalue survdiff
#' @param x Logrank test result, as returned by \code{\link[survival]{survdiff}}.
#' @param log_p Whether to return the logarithm of the p-value.
#' @param ... Ignored. Kept for S3 consistency.
#' @return p-value.
#' library(survival)
#' y <- Surv(time=1:100, event=rep(1:0, each=50))
#' groups <- rep(1:2, each=50)
#' x <- survdiff(y ~ groups)
#' 
#' # Compare p-values of implementations
#' print(x)
#' pvalue(x)
#' @seealso \code{\link{pvalue}}
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
pvalue.survdiff <- function(x, log_p=FALSE, ...){
    pchisq(x$chisq, length(x$n) - 1, lower.tail=FALSE, log.p=log_p)
}


#' Dichotomize time-to-event data
#'
#' Convert time-to-event data (typically created with the \code{\link{Surv}}
#' function) to factor or integer.
#' 
#' If no time point is given the observation times will be stripped, leaving
#' only the event types. If a time point is given
#' observations with events occurring before \code{time} will be labelled by
#' their event type,
#' observations with events occurring after \code{time} will be labelled as
#' \dQuote{no event}, and
#' observations censored before \code{time} will be considered as missing
#' information.
#' 
#' @param x \code{\link{Surv}} vector.
#' @param time Time point to dichotomize at.
#' @param to_factor Depending on the type of \code{x} the return value may be
#'   integer or factor. Set this argument to explicitly state the return type.
#' @return Integer vector or factor.
#' @seealso \code{\link{Surv}}
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
dichotomize <- function(x, time, to_factor){
    ev <- if(missing(time)){
        ifelse(is.na(x[,"time"]), NA, x[,"status"])
    } else {
        ifelse(x[,"time"] > time, 0,
               ifelse(x[,"status"] == 0, NA, x[,"status"]))
    }
    if(missing(to_factor))
        to_factor <- attr(x, "type") %in% c("mright", "mcounting")
    if(to_factor){
        labels <- Surv_event_types(x)
        factor(ev, levels=0:length(labels), labels=c("no event", labels))
    } else {
        ev
    }
}

#' Get event types of a Surv object
#'
#' @param x \code{\link{Surv}} vector.
#' @return Character vector of event types.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @noRd
Surv_event_types <- function(x){
    mstat <- attr(x, "type") %in% c("mright", "mcounting")
    if(mstat) attr(x, "states") else c("no event", "event")
}

#' Plot Surv vector
#' 
#' @method plot Surv
#' @param x \code{\link{Surv}} vector.
#' @param y Y-values.
#' @param segments Whether to draw horizontal segments.
#' @param flip Flip the plot to show time on y.
#' @param legendpos Position of legend, see \code{\link{legend}}. Set to NA or
#'   NULL to suppress legend.
#' @param ... Sent to \code{\link{plot}}.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
plot.Surv <- function(x, y, segments=TRUE, flip=FALSE, legendpos="topright", ...){
    if(missing(y)) y <- 1:length(x)

    if(flip){
        plot(y, x[,"time"], type="n", ...)
        if(segments) segments(y, 0, y, x[,"time"], col=dichotomize(x, to_factor=FALSE)+1)
        points(y, x[,"time"], pch=20, col=dichotomize(x)+1)
    } else {
        plot(x[,"time"], y, type="n", ...)
        if(segments) segments(0, y, x[,"time"], y, col=dichotomize(x, to_factor=FALSE)+1)
        points(x[,"time"], y, pch=20, col=dichotomize(x)+1)
    }
    if(!is_blank(legendpos)){
        legend(legendpos, Surv_event_types(x), pch=20,
               col=seq_along(Surv_event_types(x)))
    }
}

#' Fit Cox proportional hazards model
#' 
#' @param x Dataset.
#' @param y Response. Required if formula is missing.
#' @param formula See \code{\link{coxph}}.
#' @param ... Sent to \code{\link{coxph}}.
#' @return Fitted Cox proportional hazards model.
#' @examples
#' require(survival)
#' data(ovarian)
#' model <- fit(
#'     modeling_procedure(
#'         method = "coxph",
#'         parameter = list(formula = list(Surv(futime, fustat) ~ age))),
#'     x = ovarian, y = NULL
#' )
#' predict(model, ovarian[11:16,])
#' @seealso \code{\link{predict_coxph}}
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
fit_coxph <- function(x, y, formula = y ~ ., ...){
    nice_require("survival")
    if(any(grep("^strata(.*)$", attr(terms(formula), "term.labels"))))
        stop("Only unstratified Cox proportional hazards regression is implemented.")
    survival::coxph(formula, data=x, ...)
}

#' Predict using Cox proportional hazards model
#' 
#' @param object Fitted model, as returned by \code{\link{fit_coxph}}.
#' @param x Observations whose response is to be predicted.
#' @param ... Sent to \code{\link{predict.coxph}}.
#' @seealso \code{\link{fit_coxph}}
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
predict_coxph <- function(object, x, ...){
    nice_require("survival")
    list(prediction = predict(object, x, ...))
}


