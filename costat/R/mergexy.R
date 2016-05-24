mergexy <-
function (...) 
{
#
# Take multiple solution sets and merge them into one
#
#
# Turn our list of arguments into a list, and get its length. nl is the number
# of solution sets
#
l <- list(...)
nl <- length(l)
startpar <- endpar <- convergence <- minvar <- pvals <- NULL

#
# Process information in each solution set
#
for (i in 1:nl) {
	xy <- l[[i]]

	#
	# Check it is of right class
	#

        if (class(xy) != "csFSS") 
            stop(paste("Object: ", l[i], " is not of class csFSS and cannot be merged"))
	#
	# If it is the first one then initiate data sets
	#
        if (i == 1) {
            Ncoefs <- ncol(xy$startpar)/2
            tsx <- xy$tsx
            tsy <- xy$tsy
            tsxname <- xy$tsxname
            tsyname <- xy$tsyname
            filter.number <- xy$filter.number
            family <- xy$family
            spec.filter.number <- xy$spec.filter.number
            spec.family <- xy$spec.family
        }

        else {	# Not the first object

	    #
	    # Check new information in new object matches what was seen
	    # previously. We don't want to merge info from different types
	    # of run and different datasets.
	    #


            if (ncol(xy$startpar)/2 != Ncoefs) 
                stop(paste("Object ", l[i],
		    " is based on a different number of coefficients (", 
		    ncol(xy$startpar)/2, ") than those discovered earlier (", 
		    Ncoefs, ")"))

            if (tsxname != xy$tsxname) 
                stop(paste("Name of X series in object: ", l[i], 
		    " is different than those discovered earlier. All X names should be the same in each csFSS object"))

            if (tsyname != xy$tsyname) 
                stop(paste("Name of Y series in object: ", l[i], 
		    " is different than those discovered earlier. All Y names should be the same in each csFSS object"))

            if (any(tsx != xy$tsx)) 
                stop(paste("Entries in tsx component of object: ", 
		    l[i], " are not all the same as the one in the first object. All tsx vectors in all objects have to be the same"))

            if (any(tsy != xy$tsy)) 
                stop(paste("Entries in tsy component of object: ", 
		    l[i], " are not all the same as the one in the first object. All tsy vectors in all objects have to be the same"))

            if (filter.number != xy$filter.number) 
                stop(paste("filter.number for object: ", l[i], 
		    " is not the same as the one in the first object. All objects have to have the same filter.number"))

            if (family != xy$family) 
                stop(paste("family for object: ", l[i],
		    " is not the same as the one in the first object. All objects have to have the same family"))

            if (spec.filter.number != xy$spec.filter.number) 
                stop(paste("spec.filter.number for object: ", l[i],
		    " is not the same as the one in the first object. All objects have to have the same spec.filter.number"))

            if (spec.family != xy$spec.family) 
                stop(paste("spec.family for object: ", l[i], 
		" is not the same as the one in the first object. All objects have to have the same spec.family"))
        }
	#
	# Passed all the checks, now store the new info
	#
        startpar <- rbind(startpar, xy$startpar)
        endpar <- rbind(endpar, xy$endpar)
        convergence <- c(convergence, xy$convergence)
        minvar <- c(minvar, xy$minvar)
        pvals <- c(pvals, xy$pvals)
    }

#
# Build return object and return it
#

l <- list(startpar = startpar, endpar = endpar, convergence = convergence, 
	minvar = minvar, pvals = pvals, tsx = tsx, tsy = tsy, 
        tsxname = tsxname, tsyname = tsyname, filter.number = filter.number, 
        family = family, spec.filter.number = spec.filter.number, 
        spec.family = spec.family)
class(l) <- "csFSS"
return(l)
}
