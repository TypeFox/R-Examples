# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

setMethod("setNA",
    signature(x = "data.frame", control = "NAControl"),
    function(x, control, i = 1) {
        # initializations
        target <- getTarget(control)
        if(is.null(target)) target <- getNames(x)
        lengthTarget <- length(target)
        NArate <- getNArate(control)
        if(is(NArate, "numeric")) NArate <- NArate[i]
        else NArate <- rep(NArate[i,], length.out=lengthTarget)
        if(all(NArate == 0) || any(dim(x) == 0)) return(x)  # nothing to do
        grouping <- getGrouping(control)
        if(length(grouping) > 1) {
            stop("'grouping' must not specify more than one variable")
        }
        useGroup <- as.logical(length(grouping))  # 'grouping' supplied
        aux <- getAux(control)
#        if(length(aux) > 1) {
#            stop("'aux' must not specify more than one variable")
#        }
#        useAux <- as.logical(length(aux))  # 'aux' supplied
        if(length(aux) > 1) aux <- rep(aux, length.out=lengthTarget)
        lengthAux <- length(aux)
        useAux <- as.logical(lengthAux)  # 'aux' supplied
        intoContamination <- getIntoContamination(control)
        if(intoContamination) contaminated <- NULL
        else contaminated <- x$.contaminated
        isContaminated <- !intoContamination && !is.null(contaminated)
        # get population size and number of observations/groups to be set NA
        if(useGroup) {
            groups <- x[, grouping]  # group of each individual
            if(useAux || isContaminated) {
                Ntotal <- nrow(x)
                split <- split(1:Ntotal, getFactor(groups))
                N <- length(split)
            } else {
                uniqueGroups <- unique(groups)  # unique groups
                N <- length(uniqueGroups)  # number of groups
            }
            if(isContaminated) {
                # don't set to NA if any in the group is contaminated
                contaminated <- sapply(split, function(i) any(contaminated[i]))
            }
        } else N <- nrow(x)
        n <- ceiling(NArate * N)
        # prepare auxiliary variable, if supplied
        if(useAux) {
            auxNames <- aux
            aux <- x[, aux]
#            if(!is.numeric(aux)) {
#                stop("slot 'aux' in 'control' must specify a numeric variable.")
#            }
#            if(!all(is.finite(aux))) {
#                stop("variable definted by slot 'aux' in 'control'", 
#                    "must not contain missing or infinite values.")
#            }
#            if(any(aux < 0)) aux <- aux - min(aux)  # add small value?
            if(useGroup) {
                # use the group means (much faster than medians)
                if(lengthAux == 1) {
                    #aux <- sapply(split, function(i) median(aux[i]))
                    aux <- sapply(split, function(i) mean(aux[i]))
                } else {
                    aux <- aggregate(aux, list(getFactor(groups)), mean)[, -1]
                    names(aux) <- auxNames
                }
            }
        } else aux <- NULL
        # get indices
        if(length(n) == 1) {
            if(lengthAux > 1) {
                ind <- sapply(aux, 
                    function(a) getIndicesSetNA(N, n, a, contaminated))
            } else {
                ind <- replicate(lengthTarget, 
                    getIndicesSetNA(N, n, aux, contaminated))
            }
        } else {
            if(lengthAux > 1) {
                ind <- mapply(getIndicesSetNA, N, n, aux,
                    MoreArgs=list(contaminated=contaminated))
            } else {
                ind <- mapply(getIndicesSetNA, N, n, 
                    MoreArgs=list(aux=aux, contaminated=contaminated))
            }
        }
        if(useGroup) {
            # get indices for individuals
            if(useAux || isContaminated) {
                ind <- apply(ind, 2, 
                    function(i) {
                        ans <- logical(Ntotal)
                        ans[unlist(split[i])] <- TRUE
                        ans
                    })
            } else ind <- apply(ind, 2, function(i) groups %in% uniqueGroups[i])
        }
        ## do not append logical variables indicating the values set to NA
        ## using 'is.na' in the function for the simulation runs is much faster
        #colnames(ind) <- paste(".", target, sep="")
        # set selected values to NA and return x
        x[, target][ind] <- NA
        #existing <- intersect(names(x), colnames(ind))
        #if(length(existing)) {
        #    x[, existing] <- x[, existing] | ind[, existing]
        #    new <- setdiff(colnames(ind), existing)
        #    if(length(new)) x <- cbind(x, ind[, new, drop=FALSE])
        #    x
        #} else cbind(x, ind)
        x
    })

setMethod("setNA", 
    signature(x = "data.frame", control = "character"), 
    function(x, control, ...) {
        if(length(control) != 1) {
            stop("'control' must specify exactly one ", 
                "class extending \"VirtualNAControl\"")
        }
        if(!extends(control, "VirtualNAControl")) {
            stop(gettextf("\"%s\" does not extend \"VirtualNAControl\"", 
                    control))
        }
        setNA(x, new(control, ...))
    })

setMethod("setNA",
    signature(x = "data.frame", control = "missing"),
    function(x, control, ...) {
        setNA(x, NAControl(...))
    })


## utilities
# this is an internal function, otherwise there should be some error checking
getIndicesSetNA <- function(N, size, aux = NULL, contaminated = logical()) {
    x <- 1:N
    if(length(contaminated)) {
        nc <- !contaminated
        x <- x[nc]
        aux <- aux[nc]
    }
    ans <- logical(N)
    ans[samplex(x, size, aux)] <- TRUE
    ans
}
