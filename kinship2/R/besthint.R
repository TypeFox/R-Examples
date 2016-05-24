#
# When computer time is cheap, use this routine to get a "best" pedigree.
#  The autohint routine will rearrange sibling order, but not founder order.
# This calls autohint with every possible founder order, and finds that
#   plot with the least "stress"
#        wt[1] * number of dotted arcs
#      + wt[2] * lengths of dotted arcs
#      + wt[3] * lengths of parent-child bends
# perfect: if anything meets this stress level, keep it
# 
besthint <- function(ped, wt=c(1000, 10, 1), tolerance=0) {
    #
    # find founders married to founders
    #  the female of such pairs determines the plot order of founders
    #
    mom <- match(ped$momid, ped$id)
    dad <- match(ped$dadid, ped$id) 
    founders <- ped$id[is.na(mom) & is.na(dad)]  #founders and marry-ins
    fpair <-  !(is.na(match(ped$momid, founders)) | 
                is.na(match(ped$dadid, founders)))
    fmom  <- unique(match(ped$momid[fpair], ped$id)) #row num of founding moms

    # This function generates the permutations one after the other
    permute <- function(x) {
        n <- length(x)
        if (n==3) rbind(x[1:3], x[c(2,1,3)], x[c(3,1,2)])
        else {
            temp <- paste("cbind(x[", 1:n, "], permute(x[-", 1:n, 
                          "]))", collapse=',')
            temp <- paste("rbind(", temp, ")")
            eval(parse(text=temp))
            }
        }
    pmat <- permute(1:length(fmom))
    # Put the subsets into a random order
    #  For most pedigrees, there are several permutations that will give a
    #  tolerance or near tolerance plot.  This way we should hit one of them soon.
    pmat <- pmat[order(runif(nrow(pmat))),]

    n <- length(ped$id)
    for (perm in 1:nrow(pmat)) {
        hint <- cbind(1:n, rep(0,n))
        hint[fmom,1] <- pmat[perm,]
        
        ped$hints <- hint
        newhint <- autohint(ped)  #this fixes up marriages and such
        plist <- align.pedigree(ped, packed=TRUE, align=TRUE, 
                                width=8, hints=newhint)
        
        # Compute the error measures
        err <- rep(0.0,3)
        maxlev <- nrow(plist$nid)
        for (lev in 1:maxlev) {
            idlist <- plist$nid[lev, 1:plist$n[lev]]
            dups <- duplicated(idlist)
            if (any(dups)) {
                err[1] <- err[1] + sum(dups)
                for (i in idlist[dups]) {
                    who <- (1:length(idlist))[match(idlist, i, nomatch=0)>0]
                    err[2] <- err[2] + abs(diff(plist$pos[lev, who]))
                    }
                }
            
            # get parent-child pulls
            fam2 <- plist$fam[lev,]
            if (any(fam2>0)){
                centers <- tapply(plist$pos[lev,], fam2, mean)  #center of kids
                if (any(fam2==0)) centers <- centers[-1]
                above <- plist$pos[lev-1, sort(unique(fam2))] + .5 #parents 
                err[3] <- err[3] + sum(abs(centers-above))
                }
            }
        
        # best one so far?
        total <- sum(err * wt)
        if (perm==1 || total < besttot) {
            besttot <- total
            besthint <- newhint
            }
        if (besttot<=tolerance) break   #we needn't do better than this!
        }
    besthint
    }

                
        
