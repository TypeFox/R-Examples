`allPerms` <- function(n, control = how(), check = TRUE) {
    ## start
    v <- n
    ## expand n if a numeric or integer vector of length 1
    if((is.numeric(n) || is.integer(n)) && (length(n) == 1))
         v <- seq_len(n)
    ## number of observations in data
    n <- nobs(v)
    ## check permutation scheme and update control
    make <- getMake(control)
    if (check) {
        control2 <- control
        setMake(control2) <- FALSE
        pcheck <- check(v, control = control2, quietly = TRUE)
    }
    ## ctrl <- pcheck$control
    ## if we do copy the new updated control, we need to update to
    ## reset make
    ## ctrl <- update(ctrl, make = make)

    ## get max number of permutations
    nperms <- numPerms(v, control = control)

    ## sanity check - don't let this run away to infinity
    ## esp with type = "free"
    if(nperms > getMaxperm(control))
        stop("Number of possible permutations too large (> 'maxperm')")

    WI <- getWithin(control)
    strataP <- getStrata(control, which = "plots")
    typeW <- getType(control, which = "within")
    typeP <- getType(control, which = "plot")
    BLOCKS <- getBlocks(control)
    dimW <- getDim(control, which = "within")
    dimP <- getDim(control, which = "plots")
    mirrorW <- getMirror(control, which = "within")
    mirrorP <- getMirror(control, which = "plots")
    constantW <- getConstant(control)

    ## give a BLOCKS if non supplied - i.e. one block
    if(is.null(BLOCKS))
        BLOCKS <- factor(rep(1, n))

    ## split v by blocks
    spl <- split(seq_len(n), BLOCKS)
    nb <- length(spl) # number of blocks

    ## result object
    out <- vector(mode = "list", length = nb)

    ## null-out Blocks in control
    control2 <- control
    setBlocks(control2) <- NULL

    ## loop over blocks and return allPerms on each block
    for (i in seq_along(spl)) {
        out[[i]] <-
            doAllPerms(spl[[i]], strataP, typeW, typeP, mirrorW,
                       mirrorP, constantW, dimW, dimP, control2,
                       nperms = nperms)
    }

    ## bind all blocks together, repeating them as required
    out <- cbindAllPerms(out)

    ## bind all the blocks together
    ## out <- do.call(cbind, out) ## hmm are any of these the same shape?
    out[, unlist(spl)] <- out  ## is this being done at the doAllPerms level?

    if(!(observed <- getObserved(control))) {
        obs.v <- seq_len(n)
        obs.row <- apply(out, 1, function(x, obs.v) all(x == obs.v), obs.v)
        out <- out[!obs.row, , drop = FALSE]
        ## reduce the number of permutations to get rid of the
        ## observed ordering
        setNperm(control) <- getNperm(control) - 1
    }
    class(out) <- c("allPerms", "matrix")
    attr(out, "control") <- control
    attr(out, "observed") <- observed
    out
}


`doAllPerms` <- function(obs, strataP, typeW, typeP, mirrorW, mirrorP,
                         constantW, dimW, dimP, control, nperms) {
    n <- length(obs)

    ## subset strataP to take only the obs indices and drop the unused
    ## levels
    if (!is.null(strataP)) {
        strataP <- droplevels(strataP[obs])
    }

    ## also need to update the $strata component of control
    ## FIXME: this really should have a toplevel function to set/update
    ## sub-components of control
    ## Pl <- getPlots(control)
    ## setStrata(Pl) <- strataP
    ## setPlots(control) <- Pl
    control$plots$strata <- strataP

    ## permuting within?
    if (typeW != "none") {
        if(is.null(strataP)) {
            ## no plot-level permutations
            res <- switch(typeW,
                          free = allFree(n),
                          series = allSeries(n, nperms, mirrorW),
                          grid = allGrid(n, nperms, dimW[1],
                          dimW[2], mirrorW, constantW))
            ## use res to index original observation indices in this group
            res[] <- obs[res]
        } else {
            ## permuting within plots
            tab <- table(strataP)
            pg <- unique(tab)
            ng <-  length(tab)
            if(constantW) {
                ## same permutation in each plot
                controlW <- how(within = getWithin(control))
                nperms <- numPerms(pg, controlW)
                ord <- switch(typeW,
                              free = allFree(pg),
                              series = allSeries(pg, nperms, mirrorW),
                              grid = allGrid(pg, nperms, dimW[1],
                              dimW[2], mirrorW, constantW))
                res <- vector(mode = "list", length = ng)
                ss <- seq(0, to = prod(pg, ng-1), by = pg)
                for (i in seq_len(ng)) {
                    res[[i]] <- ord + ss[i]
                }
                ## same permutation within plots, so just cbind rather than
                ## cbindAllPerms as we don't need all combns of rows
                res <- do.call(cbind, res)
                res[] <- obs[res] ## index into the observations in this block
            } else {
                ## different permutations within plots
                nperms <- numPerms(sum(tab), control)

                if(length(pg) > 1) {
                    ## different number of observations per level of strata
                    if(typeW == "grid")
                        ## FIXME: this should not be needed once all checks are
                        ## in place in check()
                        stop("Unbalanced grid designs are not supported")
                    controlW <- how(within = getWithin(control))
                    res <- vector(mode = "list", length = ng)
                    add <- c(0, cumsum(tab)[1:(ng-1)])
                    for(j in seq_along(tab)) {
                        np <- numPerms(tab[j], controlW)
                        ord <- switch(typeW,
                                      free = allFree(tab[j]),
                                      series = allSeries(tab[j], np, mirrorW))
                        res[[j]] <- ord + add[j]
                    }
                    res <- cbindAllPerms(res)
                    res[] <- obs[res]
                } else {
                    ## same number of observations per level of strata
                    controlW <- how(within = getWithin(control))
                    np <- numPerms(pg, controlW)
                    ord <-
                        switch(typeW,
                               free = allFree(pg),
                               series = allSeries(pg, np, mirrorW),
                               grid = allGrid(pg, np, dimW[1],
                               dimW[2], mirrorW, constantW))
                    res <- vector(mode = "list", length = ng)
                    ss <- seq(0, to = prod(pg, ng-1), by = pg)
                    for(i in seq_len(ng)) {
                        res[[i]] <- ord + ss[i]
                    }
                    res <- cbindAllPerms(res)
                    res[] <- obs[res]
                }
            }
        }
    }
    ## Do we need to permute plots?
    if (!is.null(strataP) && !isTRUE(all.equal(typeP, "none"))) {
        ## permuting plots ONLY
        if(typeW == "none") {
            res <- allStrata(n, control = control)
        } else {
            ## FIXME - this need updating to work with the new code
            ## permuting blocks AND within blocks
            ## need a local CONTROL that just permutes blocks
            controlP <- how(plots = Plots(strata = strataP, type = typeP),
                            within = Within(type = "none", constant = constantW))
            ## FIXME - the above should really only need to update
            ## within as shown, not fiddle with Plots

            ## number of permutations for just the block level
            permP <- numPerms(n, control = controlP)
            ## get all permutations for the block level
            shuffP <- allStrata(n, control = controlP)
            ## copy the set of permutations for within blocks
            ## permP times - results is a list
            resP <- rep(list(res), permP)
            resP <- lapply(seq_along(resP),
                            function(k, wi, bl) {
                                t(apply(wi[[k]], 1,
                                        function(x, bl, kk) {
                                            x[bl[kk,]]
                                        }, bl = bl, kk = k))
                            },
                            wi = resP, bl = shuffP)
            res <- do.call(rbind, resP)
        }
    }
    ## some times storage.mode of res is numeric, sometimes
    ## it is integer, set to "integer" for comparisons using
    ## identical to match the observed ordering
    storage.mode(res) <- "integer"

    ## return
    res
}

## enumerate all possible permutations for a more complicated
## design
## fac <- gl(2,6)
##ctrl <- how(type = "grid", mirror = FALSE, strata = fac,
##                    constant = TRUE, nrow = 3, ncol = 2)
## ctrl <- how(strata = fac,
##                     within = Within(type = "grid", mirror = FALSE,
##                     constant = TRUE, nrow = 3, ncol = 2),
##                     blocks = Blocks(type = "free"))
## Nobs <- length(fac)
## numPerms(seq_len(Nobs), control = ctrl)
## numPerms(Nobs, control = ctrl) ## works just as well
## (tmp <- allPerms(Nobs, control = ctrl, observed = TRUE))
## (tmp2 <- allPerms(Nobs, control = ctrl))
