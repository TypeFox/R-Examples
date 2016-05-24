### ===== expert =====
###
### Internal calculations for predefined weights model
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>,
### Mathieu Pigeon <mathieu.pigeon.3@ulaval.ca>

weights <- function(nprobs, nexp, nseed, true.seed, qseed, qk, w)
{
    ## Compute set of weights.
    if (is.null(w))
        w <- rep(1/nexp, nexp)

    if (!isTRUE(all.equal(sum(w), 1)))
        stop("weights must sum to 1")

    ## Compute aggregated distribution
    breaks <- drop(matrix(qseed[, nseed + 1, 1:nexp], ncol = nexp) %*% w)

    ## Results
    structure(list(breaks = breaks, probs = qk),
              method = "Predefined Weights")
}
