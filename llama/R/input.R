input <-
function(features, performances, successes=NULL, costs=NULL, extra=NULL, minimize=T) {
    if(nrow(features) != nrow(performances) || (!is.null(successes) && (nrow(features) != nrow(successes) || nrow(performances) != nrow(successes)))) {
        warning("Different number of rows in data frames, taking only common rows.")
    }
    common = intersect(names(features), names(performances))
    if(length(common) == 0) {
        stop("Performance can't be linked to features -- no common columns!")
    }
    pnames = setdiff(names(performances), common)

    ids = factor(do.call(paste0, features[common]))
    if(length(ids) != length(levels(ids))) {
        stop("Common columns do not provide unique IDs!")
    }

    combined = merge(features, performances, by=common)
    snames = NULL
    if(!is.null(successes)) {
        commonSuccess = intersect(names(successes), common)
        successNames = setdiff(names(successes), common)
        if(length(commonSuccess) == 0) {
            stop("Successes can't be linked to rest of data -- no common columns!")
        }
        if(length(intersect(successNames, pnames)) != length(successNames)) {
            stop(paste("Successes (", paste(successNames, collapse=", "), ") and performances (", paste(pnames, collapse=", "), ") for different algorithms given!", sep=""))
        }
        names(successes) = c(commonSuccess, sapply(successNames, function(x) { paste(x, "success", sep="_") }))
        # reorder successes according to performances
        successes = successes[c(commonSuccess, sapply(successNames[order(match(successNames, pnames))], function(x) { paste(x, "success", sep="_") }))]
        combined = merge(combined, successes, by=common)
        snames = setdiff(names(successes), common)
    }

    optfun = if(minimize) { min } else { max }
    best = apply(combined, 1,
        function(x) {
            tosel = pnames
            if(!is.null(successes)) {
                nosuccs = sapply(snames[which(as.logical(x[snames]) == FALSE)], function(y) { unlist(strsplit(y, "_"))[1] })
                tosel = setdiff(pnames, nosuccs)
            }
            if(length(tosel) == 0) {
                # nothing was able to solve this instance
                return(NA)
            } else {
                perfs = as.numeric(x[tosel])
                return(tosel[which(perfs == optfun(perfs))])
            }
        })
    # simplify...
    names(best) = NULL


    retval = list(data=combined, features=setdiff(names(features), common),
        performance=pnames, success=snames, minimize=minimize,
        best=best, ids=common)
    
    if(!is.null(extra)) {
        extraNames = setdiff(names(extra), common)
        combined = merge(combined, extra, by=common)
        retval$extra = extraNames
        retval$data = combined
    }

    if(!is.null(costs)) {
        if(is.numeric(costs)) {
            # same cost for everything
            combined$cost = rep.int(costs, times=length(best))
            cnames = c("cost")
        } else if(is.data.frame(costs)) {
            commonCosts = intersect(names(costs), common)
            if(length(commonCosts) == 0) {
                stop("Costs can't be linked to rest of data -- no common columns!")
            }
            cnames = setdiff(names(costs), common)
            if(length(intersect(cnames, retval$features)) != length(retval$features)) {
                stop(paste("Costs (", paste(cnames, collapse=", "), ") and features (", paste(retval$features, collapse=", "), ") are disjoint!", sep=""))
            }
            cnames = as.vector(sapply(cnames, function(x) { paste(x, "cost", sep="_") }))
            names(costs) = c(commonCosts, cnames)
            combined = merge(combined, costs, by=common)
        } else if(is.list(costs)) {
            # cost groups
            dups = intersect(names(costs$groups), retval$features)
            if(length(dups) > 0) {
                stop(paste("Some cost groups have the same names as features:", paste(dups, collapse=", ")))
            }
            retval$costGroups = costs$groups
            commonCosts = intersect(names(costs$values), common)
            cnames = setdiff(names(costs$values), common)
            if(length(commonCosts) == 0) {
                stop("Costs can't be linked to rest of data -- no common columns!")
            }
            if(length(intersect(do.call(c, costs$groups), retval$features)) != length(retval$features)) {
                stop(paste("Cost groups (", paste(do.call(c, costs$groups), collapse=", "), ") and features (", paste(retval$features, collapse=", "), ") are disjoint!", sep=""))
            }
            if(length(intersect(names(costs$groups), cnames)) != length(names(costs$groups))) {
                stop(paste("Cost groups (", paste(do.call(c, costs$groups), collapse=", "), ") and costs (", paste(cnames, collapse=", "), ") are disjoint!", sep=""))
            }
            combined = merge(combined, costs$values, by=common)
        } else {
            stop("Invalid format for costs!")
        }
        retval$cost = cnames
        retval$data = combined
    }

    class(retval) = "llama.data"
    return(retval)
}
