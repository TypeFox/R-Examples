### ===== expert =====
###
### Internal calculations for Cooke model
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>,
### Mathieu Pigeon <mathieu.pigeon.3@ulaval.ca>

cooke <- function(nprobs, nexp, nseed, true.seed, qseed, qk, alpha)
{
    ## Compute first set of weights and calibration coefficients.
    cw <- cooke.weights(if (is.null(alpha)) 0 else alpha,
                        nprobs, nexp, nseed, true.seed,
                        qseed[, -(nseed + 1), , drop = FALSE], qk)
    w <- cw$weights

    ## If 'alpha' is NULL, the function determines confidence level
    ## that maximizes the weight given by the model to a virtual
    ## expert who would give aggregated distribution for seed
    ## variables.
    if (is.null(alpha))
    {
        ## No true optimization, here, since the function is not
        ## smooth enough. Instead, we try each value for alpha from a
        ## smart set: values just under the calibration components of
        ## each expert as obtained above.
        alpha <- c(cw$calibration - .Machine$double.eps, 1)

        ## Add virtual expert
        qseed1 <- qseed[, -(nseed + 1), , drop = FALSE]
        qseed1 <- array(c(qseed1, apply(qseed1, 2, "%*%", w)),
                        dim = c(nprobs, nseed, nexp + 1))

        ## Compute set of weights for each alpha
        w <- sapply(alpha, FUN <- function(x) cooke.weights(x, nprobs,
                                                            nexp + 1,
                                                            nseed, true.seed,
                                                            qseed1,
                                                            qk)$weights)

        wmax <- which.max(w[nexp + 1,]) # max weight for virtual expert
        alpha <- alpha[wmax]            # optimized alpha
        w <- w[-(nexp + 1), wmax]       # optimal set of weights
        w <- w / sum(w)                 # normalization
    }

    ## Compute aggregated distribution
    breaks <- drop(matrix(qseed[, nseed + 1, 1:nexp], ncol = nexp) %*% w)

    ## Results
    structure(list(breaks = breaks, probs = qk, alpha = alpha),
              method = "Cooke")
}

cooke.weights <- function(alpha, nprobs, nexp, nseed, true.seed, qseed, qk)
{
    ## Function to determine in which interquantile space fall the
    ## true values of the seed variable 'i' for each expert. Returns a
    ## (nprobs - 1) x nexp boolean matrix.
    s <- seq_len(nprobs - 1) # interquantile spaces IDs
    fun <- function(i)
        apply(matrix(qseed[, i, ], ncol = nexp), 2,
              function(v) cut(true.seed[i], v, labels = FALSE) == s)

    ## Compute the empirical distribution for each interquantile space
    ## and each expert.
    S <- array(rowSums(sapply(seq_len(nseed), fun))/nseed, c(nprobs - 1, nexp))

    ## Calibration
    calibration <- pchisq(2 * nseed * colSums(S * pmax(log(S/qk), 0)),
                  df = nseed - 1, lower.tail = FALSE)

    ## Entropy
    entropy <- colSums(apply(qseed, c(2, 3), function(x)
                         log(diff(range(x)))) +
                   colSums(qk * log(qk/apply(qseed, c(2, 3),
                                             diff)))) / nseed

    ## Weights
    nnw <- calibration * entropy * (calibration > alpha) # raw weights
    w <- nnw / sum(nnw)                 # normalized weights
    list(weights = w, calibration = calibration)
}
