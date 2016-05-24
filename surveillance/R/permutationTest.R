################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Permutation test to compare the means of paired samples
###
### Copyright (C) 2011-2012 Michaela Paul, 2013-2015 Sebastian Meyer
### $Revision: 1347 $
### $Date: 2015-05-29 11:45:51 +0200 (Fre, 29. Mai 2015) $
################################################################################

permutationTest <- function(score1, score2, nPermutation = 9999,
                            plot = FALSE, verbose = FALSE)
{
    stopifnot((nTime <- length(score1)) == length(score2),
              !is.na(score1), !is.na(score2))
    
    meanScore1 <- mean(score1)
    meanScore2 <- mean(score2)
    diffObserved <- meanScore1 - meanScore2
    
    diffMean <- replicate(nPermutation, {
        sel <- rbinom(nTime, size=1, prob=0.5)
        g1 <- (sum(score1[sel==0]) + sum(score2[sel==1]))/nTime
        g2 <- (sum(score1[sel==1]) + sum(score2[sel==0]))/nTime
        g1 - g2
    })
    
    if (isTRUE(plot)) plot <- list()
    if (is.list(plot)) {
        do.call("permtestplot", args = modifyList(
            list(permstats = diffMean, xmarks = c("observed" = diffObserved),
                 xlab = "Difference between means", ylab = "Density", main = ""),
            plot))
    }
    
    pVal <- (1+sum(abs(diffMean)>=abs(diffObserved))) / (nPermutation+1)

    pTtest <- t.test(score1, score2, paired=TRUE)$p.value
    
    if (verbose)
        cat("mean difference =", diffObserved,
            "\tp(permutation) =", pVal,
            "\tp(paired t-test) =", pTtest, "\n")
    
    list(diffObs=diffObserved, pVal.permut=pVal, pVal.t=pTtest)
}
