# Automatically generated from all.nw using noweb
makekinship <- function(famid, id, father.id, mother.id, unrelated=0) {
    n <- length(famid)
    if (length(id)    != n) stop("Mismatched lengths: famid and id")
    if (length(mother.id) != n) stop("Mismatched lengths: famid and mother.id")
    if (length(father.id) != n) stop("Mismatched lengths: famid and father.id")
    if (any(is.na(famid)))  stop("One or more subjects with missing family id")
    if (any(is.na(id)))     stop("One or more subjects with a missing id")
    if (is.numeric(famid)) {
        if (any(famid <0))      stop("Invalid family id, must be >0")
        }

    if (any(duplicated(id))) stop("Subject ids must be unique")

    famlist <- sort(unique(famid))  #same order as the counts table
    idlist <- id            # will be overwritten, but this makes it the
                            #  correct data type and length
    counts <- table(famid)
    cumcount <- cumsum(counts)    
     if (any(famid==unrelated)) {
        # Assume that those with famid of 0 are unrelated uniques
        #   (usually the marry-ins)
        temp <- match(unrelated, names(counts))
        nzero <- counts[temp]    
        counts <- counts[-temp]
        famlist <- famlist[famlist != unrelated]
        idlist[1:nzero] <- id[famid== unrelated]
        cumcount <- cumsum(counts) + nzero
        }
    else nzero <- 0
    
    mlist <- vector('list', length(counts))
    for (i in 1:length(counts)) {
        who <- (famid == famlist[i])
        if (sum(who) ==1) mlist[[i]] <- Matrix(0.5)  # family of size 1
        else {
            mlist[[i]] <- kinship(id[who], mother.id[who], father.id[who])
            }
        idlist[seq(to=cumcount[i], length=counts[i])] <- id[who]
        }

    if (nzero>0) mlist <- c(list(Diagonal(nzero)), mlist)
    kmat <- forceSymmetric(bdiag(mlist))
    dimnames(kmat) <- list(idlist, idlist)
    kmat
}
