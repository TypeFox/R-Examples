### ===== expert =====
###
### Unified interface to Cooke, Mendel-Sheridan and predefined weights
### models.
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>,
### Mathieu Pigeon <mathieu.pigeon.3@ulaval.ca>

expert <- function(x, method = c("cooke", "ms", "weights"), probs,
                   true.seed, alpha = NULL, w = NULL)
{
    nexp <- length(x)                   # nb experts
    nseed <- length(x[[1]]) - 1         # nb seed variables
    nprobs <- length(probs)             # nb quantiles

    if (any(probs < 0) || any(probs > 1))
        stop("values in 'probs' must be between 0 and 1")

    ## Extract data from the list and build an array containing
    ## data for the seed variables and the decision variable.
    qseed <- array(unlist(x), dim = c(nprobs, nseed + 1, nexp))

    ## If the 0th quantile and/or the 100th quantile were not required
    ## to the experts, they are determined by removing and adding 10 %
    ## of the smallest interval containing all quantiles given by the
    ## experts at the bounds of this interval.
    ##
    ## The original 'probs' vector is retained and given attributes to
    ## flag which of the lower and/or upper bounds is set
    ## automatically to enhance display in the summary method.
    probs.orig <- probs                 # save original values
    if (any(test <- ! c(0, 1) %in% probs))
    {
        r <- apply(qseed, 2, range)     # ranges for each variable
        lag <- diff(r) * 0.1            # lags for each variable
        lag[lag == 0] <- qseed[1, lag == 0, 1] * 0.1
        nprobs <- nprobs + sum(test)
        qseed <- apply(qseed, c(2, 3), "[", seq_len(nprobs)) # add "rows"
        if (test[1])                    # case probs = 0
        {
            attr(probs.orig, "set.lower") <- TRUE
            probs <- c(0.00, probs)
            qseed[2:nprobs, , ] <- qseed[1:(nprobs - 1), , ] # move data down
            qseed[1, , ] <- r[1, ] - lag # lower bounds
        }
        if (test[2])                    # case probs = 1
        {
            attr(probs.orig, "set.upper") <- TRUE
            probs <- c(probs, 1.00)
            qseed[nprobs, , ] <- r[2, ] + lag # upper bounds
        }
    }

    if (any(apply(qseed, c(2, 3), is.unsorted)))
        stop("'qseed' must be sorted increasingly.")

    R <- apply(qseed, 2, range)[, -(nseed + 1)]
    for (i in seq_len(nseed))
        if (true.seed[i] < R[1, i] || true.seed[i] > R[2, i])
            warning(paste("true value of seed variable ", i," out of range"))

    method <- match.arg(method)

    res <- switch(method,
                  cooke = cooke(nprobs, nexp, nseed, true.seed, qseed,
                                qk = diff(probs), alpha = alpha),
                  ms = MS(nprobs, nexp, nseed, true.seed, qseed,
                          qk = diff(probs)),
                  weights = weights(nprobs, nexp, nseed, true.seed, qseed,
                                    qk = diff(probs), w = w))

    res$nexp <- nexp
    res$nseed <- nseed
    res$quantiles <- probs.orig
    class(res) <- "expert"
    attr(res, "call") <- match.call()
    res
}
