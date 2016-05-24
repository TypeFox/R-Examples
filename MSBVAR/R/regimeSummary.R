# regimeSummary()
# Provides summaries of the transition probabilities and regime
# durations from a posterior MSBVAR object

regimeSummary <- function(x, quantiles=c(0.025,0.25,0.5,0.75,0.975))
{
    if(class(x)!="MSBVAR") {
        stop("This function can only be used with Gibbs / MCMC output from gibbs.msbvar()\n")
    }
    # Put some labels on the object so we get a sensible output
    qij <- expand.grid(1:x$h, 1:x$h)
    colnames(x$Q.sample) <- paste("q_", qij[,1], qij[,2], sep="")

    quant <- summary(x$Q.sample, quantiles=quantiles)

    # Quantiles
    cat("##############################################################\n")
    cat("Summary of MCMC draws for elements of the transition matrix Q\n")
    cat("##############################################################\n")
    print(quant)

    # Print the matrix

    cat("##############################################################\n")
    cat("Full mean transition matrix\n")
    cat("##############################################################\n")
    tmp <- matrix(quant$statistics[,1], x$h, x$h)
    colnames(tmp) <- rownames(tmp) <- paste("Regime", 1:x$h)
    print(tmp)

    # Now compute the long run regime probabilities
    N2 <- length(x$ss.sample)
    lrQ <- sapply(1:N2, function(i) {
        steady.Q(matrix(x$Q.sample[i,],x$h,x$h))})

    rownames(lrQ) <- paste("Regime", 1:x$h)

    lrQquant <- summary(mcmc(t(lrQ)), quantiles=quantiles)

    # Print
    cat("\n##############################################################\n")
    cat("Ergodic regime probabilities\n")
    cat("##############################################################\n")
    print(lrQquant)

    # Now do the regime durations.

    dur <- 1/(1-lrQ)
    rownames(dur) <- paste("Regime", 1:x$h, "duration")

    durquant <- summary(mcmc(t(dur)), quantiles=quantiles)

    cat("##############################################################\n")
    cat("Regime durations and quantiles\n")
    cat("##############################################################\n")
    print(durquant)

    invisible(list(Q.summary = quant, lrQ = lrQquant, durations=durquant))
}
