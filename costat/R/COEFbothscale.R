COEFbothscale <-
function (l, plotclustonly=FALSE, StyPval=0.05, ...)
{
#
# Check class of input object
#

if (class(l) != "csFSS")
	stop(paste("Object: ", deparse(substitute(l)), " needs to be a `csFSS' class object"))

#
# Store final coefficients in a more convenient form to use
#

endpar <- l$endpar

#
# Generate basic set of indices to reference solutions, will pare this down
# as we go through
#

ix <- 1:nrow(endpar)

#
# Do some error checking, ie can't do stuff if there are too few solns
#

if (nrow(endpar) == 1)
	stop("Can't do any comparison between solutions when there is only one solution!")

if (nrow(endpar) == 2)
	stop("No point doing scaling or hierarchical clustering with only two solutions!")

#
# Only pick solutions that had converged properly
#
endpar <- endpar[l$convergence==0,]

#
# Deal with solutions that are just the negative of others, by making all
# first coefficients have the same sign.
#
firstofeachrow <- endpar[, 1]
pon <- firstofeachrow < 0
signmult <- (-1)^pon
endpar2 <- endpar*signmult
endpar <- endpar2

#
# Extract p-values of the solutions that converged
# Choose indices and solutions that converged and were stationary
#
pvals <- l$pvals[l$convergence==0]
ix <- ix[l$convergence==0]

if (sum(pvals >= StyPval) <= 2)
	stop("There are <= 2 stationary series, not enough to proceed")

ix <- ix[pvals >= StyPval]
endpar <- endpar[pvals >=StyPval,]
dimnames(endpar) <- list(as.character(ix), NULL)

#
# Work out number of coefficients
#
lv <- ncol(endpar)/2


# Normalize rows: this is so scaling can compare properly.

mynorm <- function(x) sqrt(sum(x^2))

swst <- apply(endpar, 1, mynorm)
endpar <- sweep(endpar, 1, swst, FUN="/")

#
# Compute inter-row distances between the normalized solutions of those
# that converged and were stationary
# Then do scaling and hierarchical clustering with cmdscale
#
epd <- dist(endpar)
epscale <- cmdscale(epd)
epclust <- hclust(epd)

#
# Prepare and return answer object
#

l <- list(epscale=epscale, epclust=epclust, x=l, StyPval=StyPval)
class(l) <- "csFSSgr"
return(l)
}
