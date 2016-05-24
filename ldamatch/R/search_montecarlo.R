#' Searches for a subset satisfying the test by randomly selecting
#' subspaces of decreasing size.
#'
#' @param sspace  An ordered subject subspace: a list of vectors,
#' with one vector per group containing the corresponding subject indices.
#'
#' @inheritParams ldamatch
#'
#' @importFrom stats rbinom
search_montecarlo <- function(sspace, condition, covariates,
                              halting_test, thresh, replicates) {
    # Initializes `p` to be the minimal p.value such that the expectation
    # is that one observation will be dropped out.
    L <- vapply(sspace, length, 0)
    size <- sum(L)
    p <- 1. / size
    p.max <- 1. - p
    # Computes step for requested p-value and number of replicates.
    p.step <- (1. - 2. * p) / (replicates - 1)
    ## Computes scaling factors on p for each group.
    scaling <- L / size
    # Runs replicates.
    while (p < p.max) {
        # Initializes binary vector for included items in sample.
        is.in <- rep(TRUE, length(condition))
        # Samples to define a subset.
        for (name in names(sspace)) {
            my.group <- sspace[[name]]
            L <- length(my.group)
            n <- rbinom(1, L, p / scaling[[name]])
            if (n == 0)
                next
            is.in[my.group[1:n]] <- FALSE
        }
        # Tests and returns binary vector if anything is found.
        if (halting_test(condition[is.in], subset(covariates, is.in), thresh))
            return(is.in)
        inc(p) <- p.step
    }
    stop("Convergence failure")
}


