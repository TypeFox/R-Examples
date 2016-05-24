### ===== expert =====
###
### Internal calculations for the Mendel-Sheridan model
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>,
### Mathieu Pigeon <mathieu.pigeon.3@ulaval.ca>

MS <- function(nprobs, nexp, nseed, true.seed, qseed, qk)
{
    ## Some values and sequences occuring more than once. The second
    ## one is the column index of the quantiles given by the experts
    ## for the unknown variable.
    snexp <- seq_len(nexp)
    nseedp1 <- nseed + 1

    ## Function to determine in which joint interquantile space falls
    ## the true value of seed variable j. The 'cbind' is required in
    ## case there is only one expert.
    fun <- function(j)
        c(apply(cbind(qseed[, j, ]), 2, cut, x = true.seed[j], labels = FALSE), j)

    ## Compute calibrator for each seed variable and then the
    ## associated "sufficient" statistic.
    Z <- array(0, dim = c(rep(nprobs - 1, nexp), nseed))
    Z[t(sapply(seq_len(nseed), fun))] <- 1
    S <- apply(Z, snexp, sum)

    ## Volumes of the joint interpercentage spaces.
    a <- qk
    for (i in seq_len(nexp - 1))
        a <- outer(a, qk)

    ## Predictive distribution using the Mendel-Sheridan
    ## approximation.
    pred <- (1 + S)/(1/a + nseed)
    pred <- pred/sum(pred)

    ## Some joint interquantile spaces that were assigned a
    ## probability above are in fact impossible (no intersection
    ## between individual interquantile spaces). These must be
    ## eliminated from the output.
    ##
    ## A joint interquantile space is possible only if the largest
    ## individual space lower bound (lb) is smaller that the smallest
    ## individual space upper bound (ub).
    ##
    ## This requires to retrieve all the lower and upper bounds of the
    ## interquantile spaces for each expert. The two functions below
    ## do this for expert i.
    LOW <- function(i) cbind(qseed[, nseedp1, ])[pos[, i], i]
    HIGH <- function(i) cbind(qseed[, nseedp1,])[pos[, i] + 1, i]

    pos <- sapply(snexp, slice.index, x = a) # indexes of the spaces

    lower.bounds <- apply(sapply(snexp, LOW), 1, max)  # largest lbs
    upper.bounds <- apply(sapply(snexp, HIGH), 1, min) # smallest ubs
    is.possible <- lower.bounds < upper.bounds         # possible spaces

    pred <- pred[is.possible]           # keep only possible spaces
    pred <- pred/sum(pred)              # rescale

    ## Final results
    breaks <- c(lower.bounds[is.possible],
                tail(upper.bounds[is.possible], 1))
    structure(list(breaks = breaks, probs = pred),
              method = "Mendel-Sheridan")
}
