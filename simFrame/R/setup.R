# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

## wrapper

# name of control class can be supplied as character string
setMethod("setup",
    signature(x = "data.frame", control = "character"), 
    function(x, control, ...) {
        if(length(control) != 1) {
            stop("'control' must specify exactly one class ", 
                "inheriting from \"VirtualSampleControl\"")
        }
        if(!extends(control, "VirtualSampleControl")) {
            stop(gettextf("\"%s\" does not extend class \"VirtualSampleControl\"", 
                    control))
        }
#        setup(x, new(control, ...))
        # temporary solution: constructor for class "TwoStageControl" has 
        # arguments that aren't slots
        if(isTRUE(control == "TwoStageControl")) control <- TwoStageControl(...)
        else control <- new(control, ...)
        setup(x, control)
    })


# default control class
setMethod("setup",
    signature(x = "data.frame", control = "missing"),
    function(x, control, ...) {
        setup(x, SampleControl(...))
    })

# ---------------------------------------

## get sample setup using control class "SampleControl"

# group sampling: for methods with the argument 'x', the observation 
#                 with the first occurrence of the group will be used 
#                 as prototype (how to aggregate the data?)

setMethod("setup",
    signature(x = "data.frame", control = "SampleControl"), 
    function(x, control) {
        # initializations
        cnam <- names(x)
        design <- getCharacter(getDesign(control), cnam)
        grouping <- getCharacter(getGrouping(control), cnam)
        if(length(grouping) > 1) {
            stop("'grouping' must not specify more than one variable")
        }
        collect <- isTRUE(getCollect(control))
        groupSampling <- (length(grouping) > 0) && !collect
        collectGroups <- (length(grouping) > 0) && collect
        fun <- getFun(control)
        nam <- argNames(fun)  # argument names of 'fun'
        useX <- "x" %in% nam
        useN <- "N" %in% nam
        size <- getSize(control)
        useSize <- !is.null(size)
        prob <- getProb(control)
#        useProb <- !is.null(prob)
        haveNewProb <- is(prob, "character") || is(prob, "logical")
        if(haveNewProb) {
            prob <- getCharacter(prob, cnam)
            if(length(prob) > 1) {
                stop("'prob' must not specify more than one variable")
            }
            useProb <- length(prob) > 0
        } else useProb <- !is.null(prob) 
        dots <- getDots(control)
        k <- getK(control)
        
        # For now, the branches of the 'if' clauses are based on the 
        # corresponding previously defined methods of 'simSample'. It 
        # might be possible to increase the performance (C code?).
        
        if(length(design) > 0) {
            # -------------------
            # stratified sampling
            # -------------------
            split <- getStrataSplit(x, design)
            call <- call("mapply", FUN=fun, MoreArgs=dots, 
                SIMPLIFY=FALSE, USE.NAMES=FALSE)  # initialize call
            
            if(groupSampling) {
                # --------------
                # group sampling
                # --------------
                groupSplit <- lapply(split, function(s) x[s, grouping])
                iGroupSplit <- lapply(groupSplit, function(x) !duplicated(x))
                uniqueGroupSplit <- mapply(function(x, i) x[i], groupSplit, iGroupSplit, 
                    SIMPLIFY=FALSE, USE.NAMES=FALSE)  # unique groups in strata
                N <- sapply(uniqueGroupSplit, length)  # number of groups in strata
                if(useX) {
                    tmp <- mapply(function(x, i) x[i], split, iGroupSplit, 
                        SIMPLIFY=FALSE, USE.NAMES=FALSE)  # list of indices
                    xSplit <- lapply(tmp, function(s) x[s, , drop=FALSE])
                    call$x <- xSplit
                } else if(useN) call$N <- N
                if(useSize) {
                    if(length(size) > 1) size <- rep(size, length.out=length(N))
                    call$size <- size
                }
                if(useProb) {
                    if(haveNewProb) {
                        probSplit <- lapply(split, function(s) x[s, prob])
                        strata <- x[, design]
                        # probability weights for unique group member by strata
                        probSplit <- mapply(function(x, i) x[i], probSplit, 
                            iGroupSplit, SIMPLIFY=FALSE, USE.NAMES=FALSE)
                    } else {
                        Ngroups <- sum(N)  # number of groups
                        if(length(prob) != Ngroups) {
                            stop(gettextf("'prob' must be a vector of length %i", Ngroups))
                        }
                        # this saves time in case of multiple design variables
                        if(length(design) == 1) strata <- x[, design]
                        else strata <- unsplit(1:length(split), x[, design])
                        # if we determine unique groups in the following way, 
                        # it also works if different group members are in 
                        # different strata
                        iGroups <- unsplit(iGroupSplit, strata)
                        probSplit <- split(prob, strata[iGroups])
                        names(probSplit) <- NULL
                    }
                    call$prob <- probSplit
                    if(useSize) {
                        probSplit <- mapply(inclusionProb, probSplit, size, 
                            SIMPLIFY=FALSE, USE.NAMES=FALSE)
                    }
                    # reconstruct inclustion probabilities for individuals
                    prob <- mapply(expand, probSplit, groupSplit,           # from groups in 
                        uniqueGroupSplit, SIMPLIFY=FALSE, USE.NAMES=FALSE)  # each stratum
                    prob <- unsplit(prob, strata)  # from strata
                } else prob <- unsplit(size/N, x[, design])
                indices <- replicate(k, 
                    try(simEval(call, split, groupSplit, uniqueGroupSplit)), 
                    simplify=FALSE)  # repeated call
            } else {
                # -----------------------
                # sampling of individuals
                # groups may be collected
                # -----------------------
                N <- getStratumSizes(split, USE.NAMES=FALSE)  # stratum sizes
                if(useX) {
                    xSplit <- lapply(split, function(s) x[s, , drop=FALSE])
                    call$x <- xSplit
                } else if(useN) call$N <- N
                if(useSize) {
                    if(length(size) > 1) size <- rep(size, length.out=length(N))
                    call$size <- size
                }
                if(useProb) {
                    if(haveNewProb) prob <- x[, prob]
                    else {
                        Ntotal <- nrow(x)  # size of population
                        if(length(prob) != Ntotal) {
                            stop(gettextf("'prob' must be a vector of length %i", Ntotal))
                        }
                    }
                    probSplit <- lapply(split, function(s) prob[s])
                    call$prob <- probSplit
                    if(useSize) {
                        probSplit <- mapply(inclusionProb, probSplit, size, 
                            SIMPLIFY=FALSE, USE.NAMES=FALSE)
                        prob <- unsplit(probSplit, x[, design])
                    }
                } else prob <- unsplit(size/N, x[, design])
                # repeated call
                if(collectGroups) {
                    groups <- x[, grouping]  # group of each observation
                    groupSplit <- lapply(split, function(s) groups[s])
                    indices <- replicate(k, 
                        try(simEval(call, split, groupSplit, groupSplit)), 
                        simplify=FALSE)
                } else indices <- replicate(k, try(simEval(call, split)), simplify=FALSE)
            }
        } else {
            # -----------------
            # no stratification
            # -----------------
            call <- as.call(c(fun, dots))  # initialize call to 'fun'
            
            if(groupSampling) {
                # --------------
                # group sampling
                # --------------
                groups <- x[, grouping]  # group of each observation
                if(useX || (useProb && haveNewProb)) {
                    iGroups <- !duplicated(groups)  # logicals to extract unique groups
                    uniqueGroups <- groups[iGroups]  # unique groups
                } else uniqueGroups <- unique(groups)  # unique groups
                N <- length(uniqueGroups)  # number of groups
                if(useX) call$x <- x[iGroups, , drop=FALSE]
                else if(useN) call$N <- N
                if(useSize) {
                    if(length(size) > 1) size <- size[1]
                    call$size <- size
                }
                if(useProb) {
                    if(haveNewProb) prob <- x[iGroups, prob]
                    else if(length(prob) != N) {
                        stop(gettextf("'prob' must be a vector of length %i", N))
                    }
                    call$prob <- prob
                    if(useSize) prob <- inclusionProb(prob, size)
                    prob <- expand(prob, groups, uniqueGroups)
                } else prob <- rep.int(size/N, length(groups))
                indices <- replicate(k, 
                    try(simEval(call, groups=groups, unique=uniqueGroups)), 
                    simplify=FALSE)  # repeated call
            } else {
                # -----------------------
                # sampling of individuals
                # groups may be collected
                # -----------------------
                N <- nrow(x)  # size of population
                if(useX) call$x <- x
                else if(useN) call$N <- N
                if(useSize) {  # sample size
                    if(length(size) > 1) size <- size[1]
                    call$size <- size
                }
                if(useProb) {
                    if(haveNewProb) prob <- x[, prob]
                    else if(length(prob) != N) {
                        stop(gettextf("'prob' must be a vector of length %i", N))
                    }
                    call$prob <- prob
                    if(useSize) prob <- inclusionProb(prob, size)
                } else prob <- rep.int(size/N, N)
                # repeated call to fun
                if(collectGroups) {
                    groups <- x[, grouping]  # group of each observation
                    indices <- replicate(k, 
                        try(simEval(call, groups=groups, unique=groups)), 
                        simplify=FALSE)
                } else indices <- replicate(k, try(eval(call)), simplify=FALSE)
            }
        }
        
        # check for try errors and return 'SampleSetup' object
        # replace errors with empty index vectors: let 'runSimulation' handle 
        # the problems, this is important to debug model-based sampling
        notOK <- checkError(indices)
        if(all(notOK)) indices <- replicate(k, integer(), simplify=FALSE)
        else indices[notOK] <- integer()

        if(collectGroups) {
            # aggregate inclusion probabilities
            prob <- tapply(prob, groups, sum, simplify=FALSE)  # aggregate
            prob <- unsplit(prob, groups)  # blow up again
        }
        
        ## return 'SampleSetup' object
        SampleSetup(indices=indices, prob=prob, control=control)
    })

# ---------------------------------------

## get two-stage sample setup using control class "TwoStageControl"

setMethod("setup",
    signature(x = "data.frame", control = "TwoStageControl"), 
    function(x, control) {
        
        ## -----------
        ## first stage
        ## -----------
        
        # initializations
        cnam <- names(x)
        design <- getCharacter(getDesign(control), cnam)
        grouping <- getCharacter(getGrouping(control), cnam)
        lengthGrouping <- length(grouping)
        if(!(lengthGrouping %in% 1:2)) {
            stop("'grouping' must specify either one or two variables")
        }
        fun <- getFun(control, stage=1)
        nam <- argNames(fun)  # argument names of 'fun'
        useX <- "x" %in% nam  # use population data in call to sampling method
        useN <- "N" %in% nam  # use population size in call to sampling method
        size <- getSize(control, stage=1)
        useSize <- !is.null(size)  # use sample size in call to sampling method
        prob1 <- getProb(control, stage=1)
#        useProb <- !is.null(prob1)  # use probability weights for sampling
        haveNewProb <- is(prob1, "character") || is(prob1, "logical")
        if(haveNewProb) {
            prob1 <- getCharacter(prob1, cnam)
            if(length(prob1) > 1) {
                stop("'prob1' must not specify more than one variable")
            }
            useProb <- length(prob1) > 0
        } else useProb <- !is.null(prob1) 
        dots <- getDots(control, stage=1)
        k <- getK(control)
        
        ## computations
        PSU <- x[, grouping[1]]  # primary sampling units
        if((length(design) > 0) || useX || (useProb && haveNewProb)) {
            iPSU <- !duplicated(PSU)  # indices of unique primary sampling units
            uniquePSU <- PSU[iPSU]  # unique primary sampling units
        } else uniquePSU <- unique(PSU)
        
        if(length(design) > 0) {
            # -------------------
            # stratified sampling
            # -------------------
            
            xPSU <- x[iPSU, , drop=FALSE]  # contains values for unique PSUs
            iSplit <- getStrataSplit(xPSU, design)
            uniquePSUSplit <- lapply(iSplit, function(i) uniquePSU[i])
            N <- sapply(uniquePSUSplit, length)  # number of PSUs in strata
            # construct function call to 'fun'
            call <- call("mapply", FUN=fun, MoreArgs=dots, 
                SIMPLIFY=FALSE, USE.NAMES=FALSE)  # initialize call
            if(useX) {
                xSplit <- lapply(iSplit, function(i) xPSU[i, , drop=FALSE])
                call$x <- xSplit
            } else if(useN) call$N <- N
            if(useSize) {
                # stratification, hence 'size' should be a single integer or 
                # an integer vector with number of strata as length
                if(length(size) > 1) size <- rep(size, length.out=length(N))
                call$size <- size
            }
            # compute first stage inclusion probabilities
            if(useProb) {
                if(haveNewProb) prob1 <- xPSU[, prob1]
                else {
                    NPSU <- sum(N)  # number of PSUs
                    if(length(prob1) != NPSU) {  # check probability weights
                        stop(gettextf("'prob1' must be a vector of length %i", NPSU))
                    }
                }
                probSplit <- lapply(iSplit, function(i) prob1[i])
                names(probSplit) <- NULL
                call$prob <- probSplit  # add probability weights to function call
                if(useSize) {
                    probSplit <- mapply(inclusionProb, probSplit, size, 
                        SIMPLIFY=FALSE, USE.NAMES=FALSE)
                    # reconstruct inclustion probabilities for PSUs
                    prob1 <- unsplit(probSplit, xPSU[, design])
                }
                # reconstruct inclustion probabilities for all individuals
                prob1 <- expand(prob1, PSU, uniquePSU)
            } else {
                # first stage inclusion probabilities equal for all individuals
                # within the same stratum
                prob1 <- unsplit(size/N, x[, design])
            }
            # repeated call to sample PSUs
            # convert sampled PSUs to character strings so that the SSUs 
            # can be extracted from the list correctly with subscripts 
            indices <- replicate(k, 
                try(as.character(simEval(call, uniquePSUSplit))), 
                simplify=FALSE)
            
        } else {
            # -----------------
            # no stratification
            # -----------------
            
            N <- length(uniquePSU)  # number of PSUs
            # construct function call to 'fun'
            call <- as.call(c(fun, dots))  # initialize call
            if(useX) call$x <- x[iPSU, , drop=FALSE]
            else if(useN) call$N <- N
            if(useSize) {
                # no stratification, hence 'size' should be a single integer
                if(length(size) > 1) size <- size[1]
                call$size <- size  # add sample size to function call
            }
            # compute first stage inclusion probabilities
            if(useProb) {
                if(haveNewProb) prob1 <- x[iPSU, prob1]
                else {
                    if(length(prob1) != N) {  # check probability weights
                        stop(gettextf("'prob1' must be a vector of length %i", N))
                    }
                }
                call$prob <- prob1  # add probability weights to function call
                # compute inclusion probabilities for each individual
                if(useSize) prob1 <- inclusionProb(prob1, size)
                prob1 <- expand(prob1, PSU, uniquePSU)
            } else {
                # first stage inclusion probabilities equal for all individuals
                prob1 <- rep.int(size/N, nrow(x))
            }
            # repeated call to sample PSUs
            # convert sampled PSUs to character strings so that the SSUs 
            # can be extracted from the list correctly with subscripts 
            indices <- replicate(k, try(as.character(uniquePSU[eval(call)])), 
                simplify=FALSE)
        }
        
        # check for try errors and replace errors with empty index vectors
        notOK <- checkError(indices)
        if(all(notOK)) indices <- replicate(k, integer(), simplify=FALSE)
        else indices[notOK] <- integer()
        
        ## ------------
        ## second stage
        ## ------------
        
        # initializations
        fun <- getFun(control, stage=2)
        nam <- argNames(fun)  # argument names of 'fun' (sampling method)
        useX <- "x" %in% nam  # use population data in call to sampling method
        useN <- "N" %in% nam  # use population size in call to sampling method
        size <- getSize(control, stage=2)
        useSize <- !is.null(size)  # use sample size in call to sampling method
        prob2 <- getProb(control, stage=2)
#        useProb <- !is.null(prob2)  # use probability weights for sampling
        haveNewProb <- is(prob2, "character") || is(prob2, "logical")
        if(haveNewProb) {
            prob2 <- getCharacter(prob2, cnam)
            if(length(prob2) > 1) {
                stop("'prob2' must not specify more than one variable")
            }
            useProb <- length(prob2) > 0
        } else useProb <- !is.null(prob2) 
        dots <- getDots(control, stage=2)
        
        # get list containing the indices of SSUs in each PSU
        useSSU <- lengthGrouping == 2
        if(useSSU) {
            # secondary sampling units other than individuals (e.g., households)
            SSU <- x[, grouping[2]]  # secondary sampling units
            iSSU <- !duplicated(SSU)  # indices of unique secondary sampling units
            uniqueSSU <- SSU[iSSU]  # unique secondary sampling units
            if(useX || useProb) {
                iSplit <- split(1:length(uniqueSSU), PSU[iSSU])
                SSUbyPSU <- lapply(iSplit, function(i) uniqueSSU[i])
            } else SSUbyPSU <- split(uniqueSSU, PSU[iSSU])
        } else {
            # individuals are secondary sampling units
            SSUbyPSU <- split(1:nrow(x), PSU)
        }
        N <- sapply(SSUbyPSU, length)  # population size by PSU
        
        call <- as.call(c(fun, dots))  # initialize call to 'fun'
        if(useX) {
            # split population according to PSUs
            if(useSSU) {
                tmp <- lapply(iSplit, function(i, j) j[i], j=which(iSSU))
                xSplit <- lapply(tmp, function(i) x[i, , drop=FALSE])
            } else xSplit <- lapply(SSUbyPSU, function(i) x[i, , drop=FALSE])
        }
        if(useSize) {
            # number of observations to sample from each selected PSU
            size <- rep(size, length.out=length(N))
            names(size) <- names(SSUbyPSU)
        }
        if(useProb) {
            # second stage probability weights
            if(haveNewProb) {
                if(useSSU) prob2 <- x[iSSU, prob2]
                else prob2 <- x[, prob2]
            } else {
                NSSU <- sum(N)
                if(length(prob2) != NSSU) {  # check probability weights
                    stop(gettextf("'prob2' must be a vector of length %i", NSSU))
                }
            }
            if(useSSU) probSplit <- lapply(iSplit, function(i) prob2[i])
            else probSplit <- lapply(SSUbyPSU, function(i) prob2[i])
        }
        # function to sample SSUs within a PSU
        sampleWithinPSU <- function(i) {
            if(useX) call$x <- xSplit[[i]]
            else if(useN) call$N <- N[i]
            if(useSize) call$size <- size[i]
            if(useProb) call$prob <- probSplit[[i]]
            SSUbyPSU[[i]][eval(call)]  # sample SSUs from the currenct PSU
        }
        # function for second stage sampling from each PSU
        if(useSSU) {
            secondStage <- function(indices) {
                if(length(indices) == 0) integer()
                else try(which(SSU %in% unlist(lapply(indices, sampleWithinPSU))))
            }
        } else {
            secondStage <- function(indices) {
                if(length(indices) == 0) integer()
                else try(unlist(lapply(indices, sampleWithinPSU)))
            }
        }
        # repeated call to function for second stage sampling 
        indices <- lapply(indices, secondStage)
        # compute second stage inclusion probabilities
        if(useProb) {
            if(useSize) {
                probSplit <- mapply(inclusionProb, probSplit, size, SIMPLIFY=FALSE)
                # reconstruct inclustion probabilities
                if(useSSU) prob2 <- unsplit(probSplit, PSU[iSSU]) # for SSUs
                else prob2 <- unsplit(probSplit, PSU)  # for individuals
            }
            if(useSSU) {
                # reconstruct inclustion probabilities for all individuals
                prob2 <- expand(prob2, SSU, uniqueSSU)
            }
        } else {
            # second stage inclusion probabilities equal for all individuals
            # within the same PSU
            prob2 <- unsplit(size/N, PSU)
        }
        
        # check for try errors and replace errors with empty index vectors
        notOK <- checkError(indices)
        if(all(notOK)) indices <- replicate(k, integer(), simplify=FALSE)
        else indices[notOK] <- integer()
        
        ## return 'SampleSetup' object
        SampleSetup(indices=indices, prob=prob1*prob2, control=control)
    })

# ---------------------------------------

### utilities

# check for errors in a list (for sample setup)
checkError <- function(x) {
    # x ... list
    if(length(x)) sapply(x, function(x) class(x) == "try-error")
    else logical()
}

# function to expand a vector according to groups
# unsplit does not do the right thing here as it expects the vector to be in 
# the order of the factor levels, but we have order of first occurence
expand <- function(x, groups, unique) {
    names(x) <- as.character(unique)
    x <- x[as.character(groups)]
    unname(x)
}

# evaluate calls to sample methods with stratification and group sampling
simEval <- function(call, split, groups, unique) {
    if(!missing(split)) {
        if(missing(groups)) {
            # 'call' returns list of within-strata indices.
            tmp <- eval(call)
        } else {
            # 'groups' and 'unique' are lists
            # 'call' returns list of within-strata indices of groups. these 
            # are used to obtain the within-strata indices of individuals.
            tmp <- mapply(function(i, g, u) which(g %in% u[i]), 
                eval(call), groups, unique, SIMPLIFY=FALSE, USE.NAMES=FALSE)
        }
        # within-strata indices are turned into global indices 
        # and the resulting list is converted to a vector.
        unlist(mapply("[", split, tmp, SIMPLIFY=FALSE, USE.NAMES=FALSE))
    } else if(!missing(groups)) {  # only 'groups' is not missing
        # 'groups' and 'unique' are vectors
        # 'call' returns list of groups. these are
        # used to obtain the indices of individuals.
        which(groups %in% unique[eval(call)])
    } else eval(call)
}
