# When a random effect has two (or more) variance terms involved, then
#   we have to ensure that they can be put into the same order.  That is,
#   identical dimnames and blocksize.  
# With bdsmatrices, this isn't always so trivial, since they cannot be
#   arbitrarily reordered.  This routine does the checking, and then if
#   possible reorders them.
# "group" is a list of id's that appear in the random effects.  Any dimnames
#   that aren't in that list are not of interest to us.
# 
bdsmatrix.reconcile <- function(varlist, group) {
    ismat <- function(x) class(x) %in% c('matrix', 'bdsmatrix')
    msize <- length(group)
    #the final size of our matrices
    if(any(duplicated(group))) stop("Group index has duplicates")
    # first the easy case, only 1 matrix passed to us
    if(!is.list(varlist)) {
        if(is.function(varlist)) {
            # Someone called us with just a bdsI() call-- trivial case
            varlist <- varlist(group)
            if(!inherits(varlist, "bdsmatrix"))
                stop("Invalid function call in a varlist")
            }
        else if(ismat(varlist)) {
            #
            # Ensure that the variance matrix contain a row/col for each
            #  value in the random effect.  The 'indx' variable will
            #  contain the index for each cluster into the variance matrix.
            # Resultant coefficients are thus in the order of the dimnames
            #  of the variance matrices.  This is necessary since reorder-
            #  ing the matrix would destroy its sparse structure.
            #
            kid <- dimnames(varlist)[[1]]
            if(length(kid) == 0)
                stop("No dimnames found on a variance matrix")
            if(any(duplicated(kid)))
                stop("Duplicate dimnames in a variance matrix")
            indx <- match(group, kid)
            if(any(is.na(indx))) {
                temp <- group[is.na(indx)]
                if(length(temp) > 5)
                    temp <- temp[1:5]
                stop(paste("Group", paste(temp, collapse=' '),
                               "is in the data but not in a varlist matrix"))
                }
            # Extract the subset of varlist that corresponds to the data
            temp <- !is.na(match(kid, group))
            if(!all(temp)) {
                #some rows need to be tossed
                varlist <- varlist[temp, temp]
                }
            }
        else stop("Invalid object in a variance list")
        return(varlist)
	}

    # The interesting case -- a list was handed to us
    # First -- all the bdsmatrices must be in the same order.
    # Check for legal dimnames on all the matrices, and find out
    #  how many bdsmatrices we have.
    # Toss away any dimensions of the matrices that I don't need.
    i <- 0
    nbds <- 0
    any.matrix <- F
    # are there any ordinary matrices?
    for(j in 1:length(varlist)) {
        kmat <- varlist[[j]]
        if(ismat(kmat)) {
            i <- i + 1
            kid <- dimnames(kmat)[[1]]
            if(length(kid) == 0)
                stop("No dimnames found on a variance matrix")
            else {
                indx <- match(group, kid)
                if(any(is.na(indx))) {
                    temp <- group[is.na(indx)]
                    if(length(temp) > 5)
                        temp <- temp[1:5]
                    stop(paste("Group", paste(temp, collapse=' '),
                               "is in the data but not in a varlist matrix"))
                    }

                if(length(kid) > length(indx)) {
                    # toss unneeded rows/cols
                    indx <- sort(indx)
                    kmat <- kmat[indx, indx]
                    varlist[[j]] <- kmat
                    }

                if(inherits(kmat, "bdsmatrix")) {
                    nbds <- nbds + 1
                    blocks <- kmat@blocksize
                    rcol <- length(kmat@rmat)/nrow(kmat)
                    }
                else any.matrix <- T

                if(i == 1 || inherits(kmat, "bdsmatrix"))
                    newgroup <- kid[!is.na(match(kid, group))]
                }
            }
	}

    #
    # Now, if there are any ordinary matrices, the job is trivial
    #  Turn them all into a large bdsmatrix.  This won't happen very
    #  often, I expect.
    if(any.matrix) {
        brow <- .C("bdsmatrix_index2",
                   as.integer(1),
                   as.integer(msize),
                   rows = integer((msize * (msize + 1))/2),
                   cols = integer((msize * (msize + 1))/2))
        hash1 <- (brow$rows - 1) * msize + brow$cols
        for(i in 1:length(varlist)) {
            kmat <- varlist[[i]]
            if(is.function(kmat)) {
                # Someone called us with just a bdsI() call-- trivial case
                kmat <- kmat(group)
                if(!inherits(kmat, "bdsmatrix"))
                    stop("Invalid function call in a varlist")
                }

            kid <- dimnames(kmat)[[1]]
            indx <- match(kid, group)
 
           if(inherits(kmat, "bdsmatrix")) {
                # Turn it into a bdsmatrix with only 1 block!
                bb <- kmat@blocksize
                bsize <- sum((bb * (bb + 1))/2)
                temp <- .C("bdsmatrix_index2",
                           as.integer(length(bb)),
                           as.integer(bb),
                           rows = integer(bsize),
                           cols = integer(bsize))
                newrow <- indx[temp$rows]
                newcol <- indx[temp$cols]
                hash2 <- (pmax(newrow, newcol) - 1) * msize +
                                           pmin(newrow, newcol)
                if(length(kmat@rmat)) {
                    rdim <- dim(kmat@rmat)
                    first <- rdim[1] - rdim[2]
                    newrow <- indx[row(kmat@rmat)]
                    newcol <- indx[first + col(kmat@rmat)]
                    hash2 <- c(hash2, (pmax(newrow, newcol) -1) * msize +
                                       pmin(newrow, newcol))
                    indx <- match(hash1, hash2, nomatch = 0)
                    temp <- c(0, kmat@blocks, kmat@rmat)
                    kmat <- bdsmatrix(blocksize = msize,
                                      blocks = temp[indx + 1], 
                                      dimnames = list(group, group))
                    }
                else {
                    temp <- c(0, kmat@blocks)
                    indx <- match(hash1, hash2, nomatch = 0
                                  )
                    kmat <- bdsmatrix(blocksize = msize,
                                      blocks = temp[indx + 1], 
                                      dimnames = list(group, group))
                    }
                }
            else kmat <- bdsmatrix(blocksize = msize,
                                   blocks = c(kmat[indx, indx]),
                                   dimnames = list(group, group))
            varlist[[i]] <- kmat
            }
        return(varlist)
	}
    else group <- newgroup

    #
    # So much for the easy cases.  There exists at least one bdsmatrix, 
    #  and we need to respect it's sparseness.
    # Now, if there are 0 or 1 bdsmatrices, then group has already
    #  been reordered the way we like it.  Otherwise, we need to
    #  do the hard part -- find that bdsmatrix with the biggest blocks, 
    #  and verify that all other bdsmatrices can be coerced to fit this
    #  one's shape.
    #  
    if(nbds > 1) {
        j <- 0
        for(i in 1:length(varlist)) {
            kmat <- varlist[[i]]
            if(inherits(kmat, "bdsmatrix")) {
                if(length(kmat@rmat) > 0) {
                    # I just can't handle 2 bdsmatrices with an rmat
                    # Yes, in theory one could.
                    stop("Can't handle 2 rmats in one list")
                    }
                j <- j + 1
                kid <- dimnames(kmat)[[1]]
                indx <- match(group, kid)
                block2 <- (rep(1:length(kmat@blocksize), kmat@blocksize))
                block2 <- block2[indx]

                if(j == 1) {
                    block1 <- block2
                    save <- kid
                    # this is currently the "ruling" mat
                    blocks <- kmat@blocksize
                    }
                else {
                    ufun <- function(x) length(unique(x))
                    if(all(tapply(block2, block1, ufun) == 1)) {
                        #Every block in the prior "winner" is a strict
                        #  subset of one block in kmat.  Ergo, kmat
                        #  is larger; we have a new winner.
                        block1 <- block2
                        save <- kid
                        blocks <- kmat@blocksize
                        }
                    else if(!all(tapply(block1, block2, ufun) == 1)) {
                        # Neither is a subset of the other, which means
                        #  that the id's are the same, but the family
                        #  groupings aren't.  The user messed up.
                        stop(paste("Two variance matrices have",
                                   "incompatable structure"))
                        }
                    }
                }
            }
        group <- save
	}


    #
    # Now "group" is in the right order, and all matrices can be
    #   made to conform to it.  Make it so.
    # The "hash1" index contains the indexing for the blocks of
    #   the master matrix that we are creating.
    bsize <- sum((blocks * (blocks + 1))/2)
    brow <- .C("bdsmatrix_index2",
               as.integer(length(blocks)),
               as.integer(blocks),
               rows = integer(bsize),
               cols = integer(bsize))
    hash1 <- (brow$rows - 1) * msize + brow$cols
    for(i in 1:length(varlist)) {
        kmat <- varlist[[i]]
        if(is.function(kmat)) {
            kmat <- kmat(group)
            #create a matrix
            if(!inherits(kmat, "bdsmatrix")) 
                stop("varlist has a function that did not create a bdsmatrix")
            varlist[[i]] <- kmat
            }

        # kmat is guarranteed to be a bdsmatrix
        kid <- dimnames(kmat)[[1]]
        indx <- match(kid, group)
        if(any(indx != 1:length(indx)) || 
                    (length(kmat@blocksize) != length(blocks)) ||
                    any(kmat@blocksize != blocks)) {

            # I need to reorder it
            bb <- kmat@blocksize
            bsize <- sum((bb * (bb + 1))/2)
            temp <- .C("bdsmatrix_index2",
                       as.integer(length(bb)),
                       as.integer(bb),
                       rows = integer(bsize),
                       cols = integer(bsize))
            newrow <- indx[temp$rows]
            newcol <- indx[temp$cols]
            hash2 <- (pmax(newrow, newcol) - 1) * msize + pmin(newrow, newcol)
            indx <- match(hash1, hash2, nomatch = 0)

            if(rcol > 0) {
                if(length(kmat@rmat) > 0)
                    stop("Impossible branch! Show this message to TMT"
                         )
                # The parent we are matching has an rmat, kmat does not
                # hash3 will be the hash index for rmat
                first <- (msize - rcol)
                newrow <- rep(1:msize, rcol)
                newcol <- rep(first + 1:rcol, rep(msize, rcol))
                hash3 <- (pmax(newrow, newcol) - 1) * msize +
                               pmin(newrow, newcol)
                indx2 <- match(hash3, hash2, nomatch = 0)
                kmat <- bdsmatrix(blocksize = blocks, 
                                  blocks = c(0, kmat@blocks)[indx + 1],
                                  rmat = matrix(c(0, kmat@blocks)[indx2 + 1],
                                         ncol = rcol), 
                                  dimnames = list(group, group))
                }
            else {
                kmat@blocksize <- blocks
                kmat@blocks <- (c(0, kmat@blocks))[1 + indx]
                kmat@Dimnames <- list(group, group)
                }
            varlist[[i]] <- kmat
            }
	}
    varlist
}
