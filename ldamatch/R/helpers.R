## Helpers

`inc<-` <- function(x, value) x + value
`dec<-` <- function(x, value) x - value


#' Returns smallest halting_test() / thresh ratio, or zero if less than 1.
#'
#' @param crit  The criterion function to use, such as \code{\link{t_crit}}.
#'
#' @inheritParams ldamatch
apply_crit <- function(covariates, crit, condition, thresh) {
    min_p <- Inf
    for (i in 1:ncol(covariates)) {
        p <- crit(covariates[, i], condition)
        if (p < thresh)
            return(0.0)
        min_p <- min(p, min_p)
    }
    min_p / thresh
}


