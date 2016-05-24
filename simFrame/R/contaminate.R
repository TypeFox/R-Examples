# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

setMethod("contaminate",
    signature(x = "data.frame", control = "ContControl"),
    function(x, control, i = 1) {
        # initializations
        epsilon <- getEpsilon(control)[i]
#        if(epsilon == 0 || any(dim(x) == 0)) return(x)  # nothing to do
        if(epsilon == 0 || nrow(x) == 0) {
            x$.contaminated <- rep.int(FALSE, nrow(x))
            return(x)
        } else if(ncol(x) == 0) return(x)  # nothing to do
        target <- getTarget(control)
        if(is.null(target)) target <- getNames(x)
        grouping <- getGrouping(control)
        if(length(grouping) > 1) {
            stop("'grouping' must not specify more than one variable")
        }
        aux <- getAux(control)
        if(length(aux) > 1) {
            stop("'aux' must not specify more than one variable")
        }
        useGroup <- as.logical(length(grouping))  # 'grouping' supplied
        useAux <- as.logical(length(aux))  # 'aux' supplied
        # get population size and number of 
        # observations or groups to be contaminated
        if(useGroup) {
            groups <- x[, grouping]  # group of each individual
            if(useAux) {
                Ntotal <- nrow(x)
                split <- split(1:Ntotal, getFactor(groups))
                N <- length(split)
            } else {
                uniqueGroups <- unique(groups)  # unique groups
                N <- length(uniqueGroups)  # number of groups
            }
        } else N <- nrow(x)
        n <- ceiling(epsilon * N)
        if(useAux) {  # prepare auxiliary variable
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
                #aux <- sapply(split, function(i) median(aux[i]))
                aux <- sapply(split, function(i) mean(aux[i]))
            }
            ind <- ups(N, n, prob=aux)  # get indices (unequal prob. sampling)
        } else ind <- srs(N, n)  # get indices (simple random sampling)
        if(useGroup) {
            if(useAux) ind <- unlist(split[ind])
            else ind <- which(groups %in% uniqueGroups[ind])  # row indices
        }
        if(is(control, "DCARContControl")) {
            values <- do.call(getDistribution(control), c(n, getDots(control)))
            if(useGroup) {
                rep <- unsplit(1:n, getFactor(groups[ind]))  # replication indices
                values <- if(is.null(dim(values))) values[rep] else values[rep,]
            }
        } else if(is(control, "DARContControl")) {
            dots <- c(list(x[ind, target]), getDots(control))
            values <- do.call(getFun(control), dots)
        } else {
            stop("for user defined contamination control classes, a ", 
                "method 'contaminate(x, control, i)' needs to be defined")
        }
        values <- as.data.frame(values)
        # set contaminated values and return x
        x[ind, target] <- values
        if(is.null(x$.contaminated)) {
            contaminated <- logical(nrow(x))
            contaminated[ind] <- TRUE
            x$.contaminated <- contaminated
        } else x[ind, ".contaminated"] <- TRUE
        x
    })

setMethod("contaminate", 
    signature(x = "data.frame", control = "character"), 
    function(x, control, ...) {
        if(length(control) != 1) {
            stop("'control' must specify exactly one ", 
                "class extending \"VirtualContControl\"")
        }
        if(!extends(control, "VirtualContControl")) {
            stop(gettextf("\"%s\" does not extend \"VirtualContControl\"", 
                    control))
        }
        contaminate(x, new(control, ...))
    })

setMethod("contaminate",
    signature(x = "data.frame", control = "missing"),
    function(x, control, ...) {
        contaminate(x, DCARContControl(...))
    })
