# Automatically generated from all.nw using noweb
alignped1 <- function(x, dad, mom, level, horder, packed, spouselist){
    # Set a few constants
    maxlev <- max(level)
    lev <- level[x]
    n <- integer(maxlev)

    if (length(spouselist)==0)  spouse <- NULL
    else {
        if (any(spouselist[,1]==x)){
            sex <- 1                                  # I'm male
            sprows <- (spouselist[,1]==x & (spouselist[,4] ==spouselist[,3] |
                                            spouselist[,4] ==0))
            spouse <- spouselist[sprows, 2] #ids of the spouses
            }
        else {
            sex <- 2
            sprows <- (spouselist[,2]==x & (spouselist[,4]!=spouselist[,3] |
                                            spouselist[,4] ==0))
            spouse <- spouselist[sprows, 1]
            }
        }
    # Marriages that cross levels are plotted at the higher level (lower
    #  on the paper).
    if (length(spouse)) {
        keep <- level[spouse] <= lev
        spouse <- spouse[keep]
        sprows <- (which(sprows))[keep]
        }
    nspouse <- length(spouse)  # Almost always 0, 1 or 2
    nid <- fam <- matrix(0L, maxlev, nspouse+1)
    pos <- matrix(0.0, maxlev, nspouse +1)
    n[lev] <- nspouse +1       
    pos[lev,] <- 0:nspouse
    if (nspouse ==0) {   
        # Easy case: the "tree rooted at x" is only x itself
        nid[lev,1] <- x
        return(list(nid=nid, pos=pos, fam=fam, n=n, spouselist=spouselist))
        }
    lspouse <- spouse[spouselist[sprows,3] == 3-sex] # 1-2 or 2-1
    rspouse <- spouse[spouselist[sprows,3] == sex]   # 1-1 or 2-2
    if (any(spouselist[sprows,3] ==0)) {
        #Not yet decided spouses
        indx <- which(spouselist[sprows,3] ==0)
        nleft <- floor((length(sprows) + (sex==2))/2) #total number to left
        nleft <- nleft - length(lspouse)  #number of undecideds to the left
        if (nleft >0) {
            # JPS fixed 5/2013, don't index when nleft > length(indx)
            lspouse <- c(lspouse, spouse[indx[1:min(nleft,length(indx))]])
            indx <- indx[-(1:min(nleft,length(indx)))]
          }
        if (length(indx)) rspouse <- c(spouse[indx], rspouse)
      }

    nid[lev,] <- c(lspouse, x, rspouse)
    nid[lev, 1:nspouse] <- nid[lev, 1:nspouse] + .5  #marriages    

    spouselist <- spouselist[-sprows,, drop=FALSE]
    nokids <- TRUE   #haven't found any kids yet
    spouse <- c(lspouse, rspouse)  #reorder
    for (i in 1:nspouse) {
        ispouse <- spouse[i]
        children <- which((dad==x & mom==ispouse) | (dad==ispouse & mom==x))
        if (length(children) > 0) {
            rval1 <- alignped2(children, dad, mom, level, horder, 
                              packed, spouselist)
            spouselist <- rval1$spouselist
            # set the parentage for any kids
            #  a nuisance: it's possible to have a child appear twice, when
            #  via inbreeding two children marry --- makes the "indx" line
            #  below more complicated
            temp <- floor(rval1$nid[lev+1,])  # cut off the .5's for matching
            indx <- (1:length(temp))[match(temp,children, nomatch=0) >0]
            rval1$fam[lev+1,indx] <- i   #set the kids parentage
            if (!packed) {
                # line the kids up below the parents
                # The advantage at this point: we know that there is 
                #   nothing to the right that has to be cared for
                kidmean <- mean(rval1$pos[lev+1, indx])
                parmean <- mean(pos[lev, i + 0:1])
                if (kidmean > parmean) {
                    # kids to the right of parents: move the parents
                    indx <- i:(nspouse+1)
                    pos[lev, indx] <- pos[lev, indx] + (kidmean - parmean)
                    }
                else {
                    # move the kids and their spouses and all below
                    shift <- parmean - kidmean
                    for (j in (lev+1):maxlev) {
                        jn <- rval1$n[j]
                        if (jn>0) 
                            rval1$pos[j, 1:jn] <- rval1$pos[j, 1:jn] +shift
                        }
                    }
                }
            if (nokids) {
                rval <- rval1
                nokids <- FALSE
                }
            else {
                rval <- alignped3(rval, rval1, packed)
                }
            }
        }
    if (nokids) {
        return(list(nid=nid, pos=pos, fam=fam, n=n, spouselist=spouselist))
        }

    if (ncol(rval$nid) >= 1+nspouse) {
        # The rval list has room for me!
        rval$n[lev] <- n[lev]
        indx <- 1:(nspouse+1)
        rval$nid[lev, indx] <- nid[lev,]
        rval$pos[lev, indx] <- pos[lev,]
        }
    else {
        #my structure has room for them
        indx <- 1:ncol(rval$nid)   
        rows <- (lev+1):maxlev
        n[rows] <- rval$n[rows]
        nid[rows,indx] <- rval$nid[rows,]
        pos[rows,indx] <- rval$pos[rows,]
        fam[rows,indx] <- rval$fam[rows,]
        rval <- list(nid=nid, pos=pos, fam=fam, n=n)
        }
    rval$spouselist <- spouselist
    rval
    }
