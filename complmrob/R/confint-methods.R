#' Calculate confidence intervals
#'
#' Calculate confidence intervals for bootstrapped robust linear regression estimates with or without
#' compositional data
#'
#' @param object an object returned from \code{\link{bootcoefs}}.
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector
#'      of numbers or a vector of names. If missing, all parameters are considered.
#' @param level the confidence level required.
#' @param type the type of interval required (see the type argument of \code{\link{boot.ci}}).
#' @param ... currently ignored.
#'
#' @importFrom boot boot.ci
#' @import robustbase
#' @export
#' @describeIn confint for bootstrapped estimates of robust linear regression models for compositional data
#' @examples
#' data <- data.frame(lifeExp = state.x77[, "Life Exp"], USArrests[ , -3])
#' mUSArr <- complmrob(lifeExp ~ ., data = data)
#' bc <- bootcoefs(mUSArr, R = 200) # the number of bootstrap replicates should
#'                                  # normally be higher!
#' confint(bc, level = 0.95, type = "perc")
#'
#' ### For normal robust linear regression models ###
#' require(robustbase)
#' data(aircraft)
#'
#' mod <- lmrob(Y ~ ., data = aircraft)
#' bootEst <- bootcoefs(mod, R = 200)
#' confint(bootEst, level = 0.95, type = "perc")
confint.bccomplmrob <- function(object, parm, level = 0.95, type = c("bca", "perc", "norm", "basic", "stud"), ...) {
    type = match.arg(type);

    outtype <- switch(type,
        bca = "bca",
        perc = "percent",
        norm = "normal",
        basic = "basic",
        stud = "student",
        "percent"
    )
    ci <- lapply(object$bootres, function(boot.out) {
        boot::boot.ci(boot.out, conf = level, type = type)[[outtype]][1 , c(4, 5), drop = TRUE];
    });

    ci <- do.call(rbind, ci);
    colnames(ci) <- format.perc((1 + c(-1, 1) * level) / 2, 3);

    if(missing(parm)) {
        return(ci);
    } else {
        return(ci[parm, ]);
    }
}

#'
#' @importFrom boot boot.ci
#' @import robustbase
#' @describeIn confint for bootstrapped estimates of robust linear regression models
#' @export
confint.bclmrob <- function(object, parm, level = 0.95, type = c("bca", "perc", "norm", "basic", "stud"), ...) {
    type = match.arg(type);

    outtype <- switch(type,
        bca = "bca",
        perc = "percent",
        norm = "normal",
        basic = "basic",
        stud = "student",
        "percent"
    )

    ci <- lapply(seq_len(ncol(object$bootres$t)), function(i) {
        boot.ci(object$bootres, conf = level, type = type, index = i)[[outtype]][1, 4:5, drop = TRUE]
    });

    ci <- do.call(rbind, ci);
    colnames(ci) <- format.perc((1 + c(-1, 1) * level) / 2, 3);

    return(ci);
}

#' Simple function (just copied from the stats package) to format percentages
#'
#' @param probs the percentages
#' @param digits the number of digits
format.perc <- function (probs, digits) {
    paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits), "%");
}
