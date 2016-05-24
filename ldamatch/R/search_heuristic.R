#' At each vertex of the search graph, this takes a step which moves the
#' proportions of conditions in the subspace closer to the desired
#' (or sample) proportions.
#'
#' @param sspace  An ordered subject subspace: a list of vectors,
#' with one vector per group containing the corresponding subject indices.
#'
#' @inheritParams ldamatch
search_heuristic <- function(sspace, condition, covariates,
                             halting_test, thresh, props) {
    is.in <- rep(TRUE, length(condition))
    # Computes sample proportions to use if not specified.
    count <- table(condition)
    # Walks the search space.
    depth <- count; depth[] <- 1
    limits <- vapply(sspace, length, 0)
    while (all(depth < limits)) {
        excess <- names(count)[prop.table(count) > props]
        # If the current subspace has the same proportions of each
        # condition as the full sample, or _is_ the full sample because
        # this is the first iteration, (arbitrarily) remove the next
        # observation from the largest condition in the (sub)sample.
        if (length(excess) == 0)
            excess <- levels(condition)[which.max(count)]
        for (conditions in excess) {
            is.in[sspace[[conditions]][depth[[conditions]]]] <- FALSE
            if (halting_test(condition[is.in], subset(covariates, is.in), thresh))
                 return(is.in)
            inc(depth[[conditions]]) <- 1
            dec(count[[conditions]]) <- 1
        }
    }
    stop("Convergence failure")
}


