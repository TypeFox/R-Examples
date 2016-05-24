mgram <- function(species.d, space.d, breaks, nclass, stepsize, nperm = 1000, mrank = FALSE, nboot = 500, pboot=0.90, cboot=0.95, alternative = "two.sided", trace = FALSE)

# Mantel correlogram
# Written by Sarah C. Goslee
# 10 December 1997
# Updated 13 March 2006
#
# This function calculates a mantel correlogram for species.d (a
# lower-triangular distance matrix) based on the geographic distances 
# given in space.d (also a lower-triangular distance matrix). 
# nclass: number of distance classes 
# stepsize: width of distance classes
# nperm: number of permutations for mantel test
#
# Default is two-tailed test (H0: rM = 0; alternative = "two.sided")
# May also use one-sided test (H0: rM <= 0; alternative = "one.sided")

{
	species.d <- as.vector(species.d)
	space <- as.vector(space.d)

# use breaks if it exists.
# If nclass or stepsize aren't specified, use Sturge's rule to calculate nclass
# classes are shifted so that they don't have to start with zero
    if(missing(breaks)) {
        if (missing(nclass)) {
            if (missing(stepsize)) {
                nclass <- round(1 + 3.3 * log10(length(space.d)))
                stepsize <- (max(space.d) - min(space.d))/nclass
            }
            else {
                nclass <- round((max(space.d) - min(space.d))/stepsize)
            }
        }
        else {
            if (missing(stepsize)) {
                stepsize <- (max(space.d) - min(space.d))/nclass
            }
        }
        breaks <- seq(0, stepsize * nclass, stepsize)
    }
    else {
        nclass <- length(breaks) - 1
    }

	answer.m <- matrix(0, ncol=6, nrow=nclass)
   dimnames(answer.m) <- list(NULL, c("lag", "ngroup", "mantelr", "pval", "llim", "ulim"))
   answer.m[,4] <- rep(1, nrow(answer.m))

	for(i in 1:nclass) {
      dmin <- breaks[i]
      dmax <- breaks[i + 1]
		answer.m[i,1] <- (dmin + dmax) / 2

		space.dclass <- rep(0, length(space.d))
		space.dclass[space.d <= dmin] <- 1
		space.dclass[space.d > dmax] <- 1

		ngroup <- length(space.dclass) - sum(space.dclass)
		answer.m[i,2] <- ngroup

		if(ngroup > 0) {
			mant <- mantel(species.d ~ space.dclass, nperm=nperm, mrank=mrank, nboot=nboot, pboot=pboot, cboot=cboot)
			answer.m[i,3] <- mant[1]
			if(alternative == "two.sided")
				answer.m[i,4] <- mant[4]
			else
				answer.m[i,4] <- mant[2]
			answer.m[i,5] <- mant[5] 
			answer.m[i,6] <- mant[6]
		}			
			
		if(trace) cat(i, "\t", answer.m[i,2], "\t", answer.m[i, 3], "\n")	

	}

   results <- list(mgram = answer.m, resids = NA)
   class(results) <- "mgram"
   results
}

