# BDS test bootstrapped for BDSTEST_EWS in Early Warning Signals Author: Stephen
# R Carpenter, 22 Oct 2011 Modified by: Vasilis Dakos, January 1, 2012

#' @importFrom tseries bds.test

BDSboot <- function(X, varname, nboot, epsvec, emb) {
    # begin function
    
    StdEpsAll <- X  # name of variable for BDS
    neps <- length(epsvec)
    
    # Compute and print BDS test
    message("***********************************************", quote = FALSE)
    message(c("BDS test for ", varname), quote = FALSE)
    message(c("Embedding dimension = ", emb), quote = FALSE)
    
    BDS.data <- bds.test(StdEpsAll, m = emb, epsvec)
    
    message("BDS statistics for Nominal Data at each Epsilon", quote = FALSE)
    message(round(BDS.data$statistic, 3))
    message("P value based on standard normal", quote = FALSE)
    message(round(BDS.data$p.value, 3))
    
    # Bootstrap the BDS test
    nobs <- length(StdEpsAll)
    bootmat <- matrix(0, nrow = emb - 1, ncol = neps)  # matrix to count extreme BDS values
    for (i in 1:nboot) {
        # start bootstrap loop
        epsboot <- sample(StdEpsAll, nobs, replace = TRUE)
        BDS.boot <- bds.test(epsboot, m = emb, epsvec)
        for (im in 1:(emb - 1)) {
            # loop over embedding dimensions
            bootvec <- BDS.boot$statistic[im, ]
            N.above <- ifelse(bootvec > BDS.data$statistic[im, ], 1, 0)
            bootmat[im, ] <- bootmat[im, ] + N.above
        }
        # Report progress: if hash is removed from the next two lines, the program will
        # report each time an iteration is completed cat('iteration = ',i,' of
        # ',nboot,'\n') # flush.console()
    }  # end bootstrap loop
    
    message(" ", quote = FALSE)
    message(c("Bootstrap P estimates for ", varname), quote = FALSE)
    message(c("Bootstrap iterations = ", nboot), quote = FALSE)
    
    Pboot <- bootmat/nboot
    
    for (im in 1:(emb - 1)) {
        message(c("For embedding dimension =", im + 1), quote = FALSE)
        message(c("For epsilon = ", round(epsvec, 3), "bootstrap P = "), quote = FALSE)
        message(Pboot[im, ])
    }
    
    message("**********************************************************", quote = FALSE)
    
}  # end function 
