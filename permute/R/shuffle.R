## new version of shuffle() that allows for blocking
`shuffle` <- function(n, control = how()) {
    ## get blocking, if any
    Block <- getStrata(control, which = "blocks")
    ## If no blocking, put all samples in same block
    if(is.null(Block)) {
        Block <- factor(rep(1, n))
    } else {
        ## There was blocking so update control to remove it
        ## as we don't need it in control at the within-block
        ## permutations performed in the loop
        control <- update(control, blocks = NULL)
    }

    sn <- seq_len(n) ## sequence of samples in order of input

    ## split sn on basis of Block
    spln <- split(sn, Block)
    nb <- length(spln) ## number of blocks

    ## result list
    out <- vector(mode = "list", length = nb)

    ## loop over spln and shuffle in each split
    for(i in seq_len(nb)) {
        out[[i]] <- doShuffle(spln[[i]], control)
    }
    out <- unsplit(out, Block) ## undo the original splitting
    out
}

`doShuffle` <- function(ind, control) {
    ## collect strata at Plot level
    Pstrata <- getStrata(control, which = "plots", drop = TRUE)
    plotCTRL <- getPlots(control)
    ## ...which need to be reduced to only those for `ind`
    Pstrata <- Pstrata[ind]

    n <- length(ind)
    sn <- seq_len(n)

    ## if no strata at Plot level permute all samples using stated scheme
    if(is.null(Pstrata)) {
        perm <- shuffleNoStrata(n, control)
    } else {
        typeP <- getType(control, which = "plots")
        typeW <- getType(control, which = "within")

        ## permute Plot strata?
        if(isTRUE(all.equal(typeP, "none"))) { ## NO
            perm <- sn
        } else {                               ## YES
            flip <- runif(1L) < 0.5 ## logical, passed on & used only if mirroring
            perm <- shuffleStrata(Pstrata,
                                  type = typeP,
                                  mirror = plotCTRL$mirror,
                                  flip = flip,
                                  nrow = plotCTRL$nrow,
                                  ncol = plotCTRL$ncol)
        }

        ## permute the samples within Plot strata
        if(!isTRUE(all.equal(typeW, "none"))) { ## NOTE the `!`
            ## house keeping to track permuted strata - used later
            tab <- table(Pstrata[perm])
            levs <- names(tab) ## levels of Plot strata in this split

            ## use same permutation within each level of strata?
            withinCTRL <- getWithin(control)
            CONSTANT <- withinCTRL$constant
            if(isTRUE(CONSTANT)) {
                if(isTRUE(all.equal(typeW, "free"))) {
                    N <- unique(tab)[1L]
                    same.rand <- shuffleFree(N, N)
                } else if(isTRUE(all.equal(typeW, "series"))) {
                    start <- shuffleFree(n / length(levs), 1L)
                    flip <- runif(1L) < 0.5
                } else if(isTRUE(all.equal(typeW, "grid"))) {
                    start.row <- shuffleFree(withinCTRL$nrow, 1L)
                    start.col <- shuffleFree(withinCTRL$ncol, 1L)
                    flip <- runif(2L) < 0.5
                }
            } else {
                start <- start.row <- start.col <- flip <- NULL
            }

            ## copy perm at this stage
            tmp <- perm

            ## for each level of strata in this split, shuffle
            for(lv in levs) {
                ## must re-order strata here on basis of out as they
                ## may have been permuted above
                MATCH <- Pstrata[perm] == lv
                gr <- perm[MATCH]
                if((n.gr <- length(gr)) > 1) {
                    tmp[MATCH] <-
                        switch(typeW,
                               "free" = if(isTRUE(CONSTANT)) {
                                   gr[same.rand]
                               } else {
                                   gr[shuffleFree(n.gr, n.gr)]
                               },
                               "series" =
                               gr[shuffleSeries(seq_len(n.gr),
                                                mirror = withinCTRL$mirror,
                                                start = start, flip = flip)],
                               "grid" =
                               gr[shuffleGrid(nrow = withinCTRL$nrow,
                                              ncol = withinCTRL$ncol,
                                              mirror = withinCTRL$mirror,
                                              start.row = start.row,
                                              start.col = start.col,
                                              flip = flip)]
                               )
                }
            }
            perm <- tmp
        }
    }
    ind[perm]
}
