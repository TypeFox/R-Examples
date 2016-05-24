## Halting tests.

#' A univariate halting test using the t-test (t.test).
#'
#' @inheritParams ldamatch
#'
#' @export
t_halt <- function(condition, covariates, thresh) {
    apply_crit(covariates, t_crit, condition, thresh)
}


#' Criterion function for ad_halt.
#'
#' @param covariate     A vector containing a covariate to match the conditions on.
#'
#' @inheritParams ldamatch
t_crit <- function(covariate, condition) {
    stats::t.test(covariate ~ condition)$p.value
}


#' A univariate halting test using the Wilcoxon test (wilcox.test).
#'
#' @inheritParams ldamatch
#'
#' @export
U_halt <- function(condition, covariates, thresh) {
    apply_crit(covariates, U_crit, condition, thresh)
}


#' Criterion function for U_halt.
#'
#' @inheritParams t_crit
U_crit <- function(covariate, condition) {
    stats::wilcox.test(covariate ~ condition)$p.value
}


#' A univariate halting test using the Anderson-Darling test (kSamples::ad.test).
#'
#' @inheritParams ldamatch
#'
#' @export
ad_halt <- function(condition, covariates, thresh) {
    apply_crit(covariates, ad_crit, condition, thresh)
}


#' Criterion function for ad_halt.
#'
#' @inheritParams t_crit
#'
#' @import kSamples
ad_crit <- function(covariate, condition) {
    lkS <- kSamples::ad.test(
        split(covariate, condition), method = get("AD_METHOD", .ldamatch_globals),
        Nsim = get("AD_NSIM", .ldamatch_globals))
    lkS$ad[[get("AD_VERSION", .ldamatch_globals), 3]]
}


#' A univariate halting test using the Kolmogorov-Smirnov Test (ks.test).
#'
#' Note that the condition must have two levels.
#' The null hypothesis is that they are drawn from the same distribution,
#' so we want this to be *likely*. It is up to you to choose an
#' appropriate threshold.
#'
#' @inheritParams ldamatch
#'
#' @import RUnit
#'
#' @export
ks_halt <- function(condition, covariates, thresh) {
    RUnit::checkEquals(nlevels(condition), 2, "condition must have two levels")
    logical_condition <- (condition == levels(condition)[[1]])
    apply_crit(covariates, ks_crit, logical_condition, thresh)
}

#' Criterion function for ks_halt.
#'
#' @param logical_condition  A logical vector separating the groups.
#'
#' @inheritParams t_crit
ks_crit <- function(covariate, logical_condition) {
    stats::ks.test(covariate[logical_condition],
                   covariate[!logical_condition])$p.value
}


#' A multivariate halting test (appropriate when nlevels(condition) > 2)).
#'
#' @inheritParams ldamatch
#'
#' @importFrom stats manova
#'
#' @export
wilks_halt <- function(condition, covariates, thresh) {
    p <- min(summary(manova(covariates ~ condition), test = "Wilks")$stats[1, 6])
    if (p < thresh)
        return(0.0)
    p / thresh
}


#' Creates halt test from multiple tests.
#'
#' The created halt test function returns the smallest p-value / threshold
#' ratio of the values produced by the supplied halt test functions.
#'
#' @param halting_tests   A vector of halt test functions with the signature:
#'                        halting_test(condition, covariates, thresh).
#'                        For the meaning of the parameters see
#'                        \code{\link{ldamatch}}.
#' @return A function that returns the minimum of all halting test values;
#'         the threshold value supplied to it is recycled for the individual
#'         functions.
#'
#' @import RUnit
#'
#' @export
create_halting_test <- function(halting_tests) {
    RUnit::checkTrue(length(halting_tests) > 1 && all(vapply(
        halting_tests, function(fn) length(formals(fn)) == 3, TRUE)),
        "halting_tests must receive a vector of functions with 3 parameters")

    function(condition, covariates, threshes) {
        min_ratio <- Inf
        for (i in seq_along(halting_tests)) {
            ratio <- halting_tests[[i]](condition, covariates, threshes[[i]])
            if (!ratio)
                return(0.0)
            min_ratio <- min(ratio, min_ratio)
        }
        min_ratio
    }
}
