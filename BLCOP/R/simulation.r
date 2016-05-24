###############################################################################
# Mango Solutions, Chippenham SN14 0SQ 2008
# sampleFrom
# Author: Francisco
###############################################################################
# DESCRIPTION: Function to sample from a distribution or mvdistribution class object
# KEYWORDS: utilities
###############################################################################

sampleFrom <- function
(
    dstn,               #  distribution or mvdistribution class object
    n = 1               #  number of samples to generate
)
{
    .assertClass(dstn, c("distribution", "mvdistribution"))
    sampleFun <- match.fun(paste("r", dstn@RName, sep = ""))
    do.call(sampleFun, c(n, as.list(dstn@parameters)))
}


# empirical CDF  utility function

.empCDF <- function(x, ordered = FALSE)
{
    if(!ordered)
        x <- sort(x, decreasing = FALSE)
    probs <- seq(from = 0, to = 1, along = x)
    approxfun(x, probs)
}

# empirical quantile utility function

.empQuantile <- function(x, ordered = FALSE)
{
    if(!ordered)
        x <- sort(x, decreasing = FALSE)
    probs <- seq(from = 0, to = 1, along = x)
    approxfun(probs, x)
}                                                     