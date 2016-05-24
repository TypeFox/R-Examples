# Automatically generated from all.nw using noweb
kindepth <- function(id, dad.id, mom.id, align=FALSE) {
    if (class(id)=='pedigree' || class(id)=='pedigreeList') {
        didx <- id$findex
        midx <- id$mindex
        n <- length(didx)
        } 
    else {
        n <- length(id)
        if (missing(dad.id) || length(dad.id) !=n)
            stop("Invalid father id")
        if (missing(mom.id) || length(mom.id) !=n)
            stop("Invalid mother id")
        midx <- match(mom.id, id, nomatch=0) # row number of my mom
        didx <- match(dad.id, id, nomatch=0) # row number of my dad
        }
    if (n==1) return (0)  # special case of a single subject 
    parents <- which(midx==0 & didx==0)  #founders

    depth <- rep(0,n)
    # At each iteration below, all children of the current "parents" are
    #    labeled with depth 'i', and become the parents of the next iteration
    for (i in 1:n) {
        child  <- match(midx, parents, nomatch=0) +
                  match(didx, parents, nomatch=0)

        if (all(child==0)) break
        if (i==n) 
            stop("Impossible pedegree: someone is their own ancestor")

        parents <- which(child>0) #next generation of parents
        depth[parents] <- i
        }
    if (!align) return(depth)

    chaseup <- function(x, midx, didx) {
        new <- c(midx[x], didx[x])  # mother and father
        new <- new[new>0]
        while (length(new) >1) {
            x <- unique(c(x, new))
            new <- c(midx[new], didx[new])
            new <- new[new>0]
            }
        x
        }
            
    dads <- didx[midx>0 & didx>0]   # the father side of all spouse pairs
    moms <- midx[midx>0 & didx>0]
    # Get rid of duplicate pairs
    dups <- duplicated(dads + moms*n)
    if (any(dups)) {
        dads <- dads[!dups]
        moms <- moms[!dups]
        }
    npair<- length(dads)
    done <- rep(FALSE, npair)  #couples that are taken care of
    while (TRUE) {
        pairs.to.fix <- (1:npair)[(depth[dads] != depth[moms]) & !done]
        if (length(pairs.to.fix) ==0) break
        temp <- pmax(depth[dads], depth[moms])[pairs.to.fix]
        who <- min(pairs.to.fix[temp==min(temp)])  # the chosen couple
        
        good <- moms[who]; bad <- dads[who]
        if (depth[dads[who]] > depth[moms[who]]) {
            good <- dads[who]; bad <- moms[who]
            }
        abad  <- chaseup(bad,  midx, didx)
        if (length(abad) ==1 && sum(c(dads,moms)==bad)==1) {
            # simple case, a solitary marry-in
            depth[bad] <- depth[good]
            }
        else {
            agood <- chaseup(good, midx, didx)  #ancestors of the "good" side
            # For spouse chasing, I need to exclude the given pair
            tdad <- dads[-who]
            tmom <- moms[-who]
            while (1) {
                # spouses of any on agood list
                spouse <- c(tmom[!is.na(match(tdad, agood))],
                            tdad[!is.na(match(tmom, agood))])
                temp <- unique(c(agood, spouse))
                temp <- unique(chaseup(temp, midx, didx)) #parents
                kids <- (!is.na(match(midx, temp)) | !is.na(match(didx, temp)))
                temp <- unique(c(temp, (1:n)[kids & depth <= depth[good]]))
                if (length(temp) == length(agood)) break
                else agood <- temp
                }

            if (all(match(abad, agood, nomatch=0) ==0)) {
                # shift it down
                depth[abad] <- depth[abad] + (depth[good] - depth[bad])
                #
                # Siblings may have had children: make sure all kids are
                #   below their parents.  It's easiest to run through the
                #   whole tree
                for (i in 0:n) {
                    parents <- which(depth==i)
                    child <- match(midx, parents, nomatch=0) +
                             match(didx, parents, nomatch=0)
                    if (all(child==0)) break
                    depth[child>0] <- pmax(i+1, depth[child>0])
                    }
                }
            }
        done[who] <- TRUE
        }
    if (all(depth>0)) stop("You found a bug in kindepth's alignment code!")
    depth
    }
