#
# Repack IBD data as a Matrix object
#
ibdMatrix <- function(id1, id2, x, idmap, diagonal) {
        
    if (!is.null(ncol(id1)) && ncol(id1)==1) id1 <- id1[,1]
    if (!is.null(ncol(id1))) {
        # can be a matrix or a data frame
        if (ncol(id1)!=3)
            stop("Argument id1 is a matrix or dataframe  but does not have 3 columns")
        if (!missing(id2))
            stop("First argument is a matrix or dataframe, but id2 argument is present")
        if (!missing(x))
            stop("First argument is a matrix or dataframe, but x argument is present")
        id2 <- id1[,2]
        x   <- id1[,3]
        id1 <- id1[,1]
    }

    # Reset each pair such that id1 <= id2
    temp <- pmin(id1, id2)
    id2  <- pmax(id1, id2)
    id1  <- temp

    # Add diagonal elements
    if (!missing(diagonal) && diagonal !=0) {
        idlist <- unique(c(id1, id2))
        id1 <- c(id1, idlist)
        id2 <- c(id2, idlist)
        x <- c(x, rep(diagonal, length(idlist)))
    }

    # Toss away any zeros and duplicates
    keep <- (x != 0 & !duplicated(cbind(id1, id2)))
    if (!all(keep)) {
        id1 <- id1[keep]
        id2 <- id2[keep]
        x   <- x[keep]
        }


    # If the set of ids is a list of integers 1-n for some n, we'll assume
    #  the the data is already in the optimal order.  Otherwise figure
    #  out families.  I expect the latter to happen rarely to never.
    temp <- sort(unique(id1, id2))
    maxid <- max(id1, id2)
    if (maxid != as.integer(maxid) || length(temp) != maxid ||
        any(is.na(match(temp, 1:maxid)))) {
        # drat, need to figure out family blocks
        idlist <- sort(unique(c(id1,id2)))
        nid <- length(idlist)
        id1 <- match(id1, idlist)
        id2 <- match(id2, idlist)
        famid <- 1:nid  #everyone a singleton
        indx1 <- sort(unique(id2))  #the result vector for remap1 below
        indx2 <- sort(unique(id1))
        while (1) {
            remap1 <- tapply(famid[id1], id2, min)  #map each id2 to min id1
            if (all(famid[indx1] == remap1)) break
            famid[indx1] <- remap1
            remap2 <- tapply(famid[id2], id1, min)
            famid[indx2] <- remap2
        }
        remap <- (1:nid)[order(famid)]  #reordering of subjects
        id1 <- match(id1, remap)
        id2 <- match(id2, remap)
        idlist <- idlist[remap]
    }
    else idlist <- 1:maxid
    
    # dimid will be the dimnames
    if (missing(idmap)) dimid <- idlist
    else {
        if (!is.null(dim(idmap)) || ncol(idmap) !=2)
            stop("idmap must have 2 columns")
        temp <- match(idlist, idmap[,1])
        if (any(is.na(temp)))
            stop("Values appear in id1 or id2 that are not in idmap")
        dimid <- idmap[temp,2]
        }

    sparseMatrix(i=id1, j=id2, x=x, symmetric=TRUE, 
                 dimnames=list(dimid, dimid))
}
