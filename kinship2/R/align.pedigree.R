# Automatically generated from all.nw using noweb
align.pedigree <- function(ped, packed=TRUE, width=10,
                           align=TRUE, hints=ped$hints) {
    if (class(ped)== 'pedigreeList') {
        nped <- length(unique(ped$famid))
        alignment <- vector('list', nped)
        for (i in 1:nped) {
            temp <- align.pedigree(ped[i], packed, width, align)
            alignment[[i]] <- temp$alignment
            }
        ped$alignment <- alignment
        class(ped) <- 'pedigreeListAligned'
        return(ped)
        }
    
    if (is.null(hints)) {
      hints <- try({autohint(ped)}, silent=TRUE)
      if(class(hints)=="try-error") hints <- list(order=1:dim(ped))
    } else {
      hints <- check.hint(hints, ped$sex)
    }
    
    n <- length(ped$id)
    dad <- ped$findex; mom <- ped$mindex  #save typing
    if (any(dad==0 & mom>0) || any(dad>0 & mom==0))
            stop("Everyone must have 0 parents or 2 parents, not just one")
    level <- 1 + kindepth(ped, align=TRUE)

    horder <- hints$order   # relative order of siblings within a family

    if (is.null(ped$relation)) relation <- NULL
    else  relation <- cbind(as.matrix(ped$relation[,1:2]), 
                            as.numeric(ped$relation[,3]))

    if (!is.null(hints$spouse)) { # start with the hints list
        tsex <- ped$sex[hints$spouse[,1]]  #sex of the left member
        spouselist <- cbind(0,0,  1+ (tsex!='male'), 
                            hints$spouse[,3])
        spouselist[,1] <- ifelse(tsex=='male', hints$spouse[,1], hints$spouse[,2])
        spouselist[,2] <- ifelse(tsex=='male', hints$spouse[,2], hints$spouse[,1])
        }
    else spouselist <- matrix(0L, nrow=0, ncol=4)

    if (!is.null(relation) && any(relation[,3]==4)) {
        # Add spouses from the relationship matrix
        trel <- relation[relation[,3]==4,,drop=F]
        tsex <- ped$sex[trel[,1]]
        trel[tsex!='male',1:2] <- trel[tsex!='male',2:1]
        spouselist <- rbind(spouselist, cbind(trel[,1],
                                              trel[,2],
                                              0,0))
        }
    if (any(dad>0 & mom>0) ) {
        # add parents
        who <- which(dad>0 & mom>0)
        spouselist <- rbind(spouselist, cbind(dad[who], mom[who], 0, 0))
        }

    hash <- spouselist[,1]*n + spouselist[,2]
    spouselist <- spouselist[!duplicated(hash),, drop=F]
    noparents <- (dad[spouselist[,1]]==0 & dad[spouselist[,2]]==0)
     ##Take duplicated mothers and fathers, then founder mothers
    dupmom <- spouselist[noparents,2][duplicated(spouselist[noparents,2])] #Founding mothers with multiple marriages
    dupdad <- spouselist[noparents,1][duplicated(spouselist[noparents,1])] #Founding fathers with multiple marriages
    foundmom <- spouselist[noparents&!(spouselist[,1] %in% c(dupmom,dupdad)),2] # founding mothers
    founders <-  unique(c(dupmom, dupdad, foundmom))    
    founders <-  founders[order(horder[founders])]  #use the hints to order them
    rval <- alignped1(founders[1], dad, mom, level, horder, 
                              packed=packed, spouselist=spouselist)

    if (length(founders)>1) {
        spouselist <- rval$spouselist
        for (i in 2:length(founders)) {
            rval2 <- alignped1(founders[i], dad, mom,
                               level, horder, packed, spouselist)
            spouselist <- rval2$spouselist
            rval <- alignped3(rval, rval2, packed)
            }
        }
    #
    # Unhash out the spouse and nid arrays
    #
    nid    <- matrix(as.integer(floor(rval$nid)), nrow=nrow(rval$nid))
    spouse <- 1L*(rval$nid != nid)
    maxdepth <- nrow(nid)

    # For each spouse pair, find out if it should be connected with
    #  a double line.  This is the case if they have a common ancestor
    ancestor <- function(me, momid, dadid) {
        alist <- me
        repeat {
            newlist <- c(alist, momid[alist], dadid[alist])
            newlist <- sort(unique(newlist[newlist>0]))
            if (length(newlist)==length(alist)) break
            alist <- newlist
            }
        alist[alist!=me]
        }
    for (i in (1:length(spouse))[spouse>0]) {
        a1 <- ancestor(nid[i], mom, dad)
        a2 <- ancestor(nid[i+maxdepth],mom, dad)  #matrices are in column order
        if (any(duplicated(c(a1, a2)))) spouse[i] <- 2
        }
    if (!is.null(relation) && any(relation[,3] < 4)) {
        twins <- 0* nid
        who  <- (relation[,3] <4)
        ltwin <- relation[who,1]
        rtwin <- relation[who,2]
        ttype <- relation[who,3]
        
        # find where each of them is plotted (any twin only appears
        #   once with a family id, i.e., under their parents)
        ntemp <- ifelse(rval$fam>0, nid,0) # matix of connected-to-parent ids
        ltemp <- (1:length(ntemp))[match(ltwin, ntemp, nomatch=0)]
        rtemp <- (1:length(ntemp))[match(rtwin, ntemp, nomatch=0)]
        twins[pmin(ltemp, rtemp)] <- ttype
        }
    else twins <- NULL
    if ((is.numeric(align) || align) && max(level) >1) 
        pos <- alignped4(rval, spouse>0, level, width, align)
    else pos <- rval$pos

    if (is.null(twins))
         list(n=rval$n, nid=nid, pos=pos, fam=rval$fam, spouse=spouse)
    else list(n=rval$n, nid=nid, pos=pos, fam=rval$fam, spouse=spouse, 
                  twins=twins)
    }
