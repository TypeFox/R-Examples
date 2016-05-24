# Automatically generated from all.nw using noweb
kinship <- function(id, ...) {
    UseMethod('kinship')
    }

kinship.default <- function(id, dadid, momid, sex, chrtype="autosome", ...) {
    chrtype <- match.arg(casefold(chrtype), c("autosome", "x"))
    if (any(duplicated(id))) stop("All id values must be unique")
    n <- length(id)
    pdepth <- kindepth(id, dadid, momid)
    if (chrtype == "autosome") {
        if (n==1) 
            return(matrix(.5,1,1, dimnames=list(id, id)))

        kmat <- diag(c(rep(.5, n), 0))  #founders

        mrow <- match(momid, id, nomatch=n+1) #row number of the mother
        drow <- match(dadid, id, nomatch=n+1) #row number of the dad 
        ## When all unrelateds, pdepth all=0. 
        ## Put c(1,) to make guard from iter 1:0
        for (depth in 1:max(c(1,pdepth))) {
            for (j in  (1:n)[pdepth==depth]) {
                kmat[,j] <-kmat[j,] <- (kmat[mrow[j],]  + kmat[drow[j],]) /2
                kmat[j,j] <- (1 + kmat[mrow[j], drow[j]]) /2
            }
        }
    }
    else if (chrtype == "x") {
        if (missing(sex) || length(sex) !=n) 
            stop("invalid sex vector")
        #1 = female, 2=male
        if (n==1) 
            return(matrix(ifelse(sex>2,sex/2,NA), 1,1, dimnames=list(id, id)))

        # kmat <- diag(c((3-sex)/2, 0)) #founders
        kmat <- diag(ifelse(sex>2, NA, c((3-sex)/2, 0)))
        mrow <- match(momid, id, nomatch=n+1) #row number of the mother
        drow <- match(dadid, id, nomatch=n+1) #row number of the dad 

        for (depth in 1:max(c(1,pdepth))) {
            for (j in (1:n)[pdepth==depth]) {
                if (sex[j] ==1) {
                    kmat[,j] <- kmat[j,] <- kmat[mrow[j],]
                    kmat[j,j]<- 1
                } 
                else if(sex[j] == 2) {
                    kmat[,j] <-kmat[j,] <- (kmat[mrow[j],]  + kmat[drow[j],]) /2
                    kmat[j,j] <- (1 + kmat[mrow[j], drow[j]]) /2
                } 
                else {
                    kmat[,j] <-kmat[j,] <- NA
                    kmat[j,j] <- NA 
                }
            }
        }
    }
    kmat <- kmat[1:n,1:n]
    dimnames(kmat) <- list(id, id)
    kmat
}
kinship.pedigree <- function(id, chrtype="autosome", ...) {
    chrtype <- match.arg(casefold(chrtype), c("autosome", "x"))
    if (any(duplicated(id$id))) stop("All id values must be unique")
    n <- length(id$id)
    pdepth <- kindepth(id)
    
    # Are there any MZ twins to worry about?
    havemz <- FALSE
    if (!is.null(id$relation) && any(id$relation$code=="MZ twin")) {
        havemz <- TRUE
        temp <- which(id$relation$code=="MZ twin")
        ## drop=FALSE added in case only one MZ twin set
        mzmat <- as.matrix(id$relation[,c("indx1", "indx2")])[temp,,drop=FALSE]
        mzgrp <- 1:max(mzmat) #everyone starts in their own group
        # The loop below will take k-1 iterations for a set labeled as
        #   (k-1):k, ..., 4:3, 3:2, 2:1;  this is the worst case.
        while(1) {
            if (all(mzgrp[mzmat[,1]] == mzgrp[mzmat[,2]])) break
            for (i in 1:nrow(mzmat)) 
                mzgrp[mzmat[i,1]] <- mzgrp[mzmat[i,2]] <- min(mzgrp[mzmat[i,]])
            }
        mzindex <- cbind(unlist(tapply(mzmat, mzgrp[mzmat], function(x) {
                                        z <- unique(x)
                                        rep(z, length(z))})),
                         unlist(tapply(mzmat, mzgrp[mzmat], function(x) {
                                        z <- unique(x)
                                        rep(z, each=length(z))})))
        mzindex <- mzindex[mzindex[,1] != mzindex[,2],]
    }
    
    if (chrtype == "autosome") {
        if (n==1) 
            return(matrix(.5,1,1, dimnames=list(id$id, id$id)))

        kmat <- diag(c(rep(.5, n), 0))  #founders
        mrow <- ifelse(id$mindex ==0, n+1, id$mindex)
        drow <- ifelse(id$findex ==0, n+1, id$findex)

        for (depth in 1:max(pdepth)) {
            indx <- which(pdepth == depth)
            kmat[indx,] <- (kmat[mrow[indx],] + kmat[drow[indx],]) /2
            kmat[,indx] <- (kmat[,mrow[indx]] + kmat[,drow[indx]]) /2
            for (j in indx) kmat[j,j] <- (1 + kmat[mrow[j], drow[j]])/2
            if (havemz) kmat[mzindex] <- (diag(kmat))[mzindex[,1]]
        }
    }
    else if (chrtype == "x") {
        sex <- as.numeric(id$sex) # 1 = female, 2=male
        if (n==1) 
            return(matrix(sex/2, 1,1, dimnames=list(id$id, id$id)))

        kmat <- diag(c((3-sex)/2, 0))  #1 for males, 1/2 for females
        mrow <- ifelse(id$mindex ==0, n+1, id$mindex)
        drow <- ifelse(id$findex ==0, n+1, id$findex)

        for (depth in 1:max(pdepth)) {
            for (j in (1:n)[pdepth==depth]) {
                if (sex[j] ==1) {
                    kmat[,j] <- kmat[j,] <- kmat[mrow[j],]
                    kmat[j,j]<- 1
                }
                else if(sex[j]==2) {
                    kmat[,j] <-kmat[j,] <- (kmat[mrow[j],]  + kmat[drow[j],]) /2
                    kmat[j,j] <- (1 + kmat[drow[j],mrow[j]]) /2
                } else {
                  kmat[,j] <-kmat[j,] <- NA
                   kmat[j,j] <- NA
                }
            if (havemz) kmat[mzindex] <- (diag(kmat))[mzindex[,1]]
            }
        }
    }
    kmat <- kmat[1:n,1:n]
    dimnames(kmat) <- list(id$id, id$id)
    kmat
}
kinship.pedigreeList <- function(id, chrtype="autosome", ...) {
    famlist <- unique(id$famid)
    nfam <- length(famlist)
    matlist <- vector("list", nfam)
    idlist  <- vector("list", nfam) #the possibly reorderd list of id values
   
    for (i in 1:length(famlist)) {
        tped <- id[i]  #pedigree for this family
        temp <- try(kinship(tped, chrtype=chrtype, ...), silent=TRUE)
        if (class(temp)=="try-error") 
            stop(paste("In family", famlist[i], ":", temp))
        else matlist[[i]] <- as(forceSymmetric(temp), "dsCMatrix")
        idlist[[i]] <- tped$id
    }

    result <- bdiag(matlist)
    if (any(duplicated(id$id))) 
        temp <-paste(rep(famlist, sapply(idlist, length)),
                     unlist(idlist), sep='/') 
    else temp <- unlist(idlist)
        
    dimnames(result) <- list(temp, temp)
    result
}
oldkinship <- function(id, ...) {
    UseMethod('oldkinship')
    }

oldkinship.default <- function(id, dadid, momid, ...) {
    n <- length(id)
    if (n==1) 
        return(matrix(.5,1,1, dimnames=list(id, id)))
    if (any(duplicated(id))) stop("All id values must be unique")
    kmat <- diag(n+1) /2
    kmat[n+1,n+1]    <- 0 

    pdepth <- kindepth(id, dadid, momid)
    mrow <- match(momid, id, nomatch=n+1) #row number of the mother
    drow <- match(dadid, id, nomatch=n+1) #row number of the dad 

    for (depth in 1:max(pdepth)) {
        indx <- (1:n)[pdepth==depth]
        for (i in indx) {
            mom <- mrow[i]
            dad <- drow[i]
            kmat[i,]  <- kmat[,i] <- (kmat[mom,] + kmat[dad,])/2
            kmat[i,i] <- (1+ kmat[mom,dad])/2
            }
        }
    
    kmat <- kmat[1:n,1:n]
    dimnames(kmat) <- list(id, id)
    kmat
    }

oldkinship.pedigree <- function(id, ...) {
    n <- length(id$id)
    if (n==1) 
        return(matrix(.5,1,1, dimnames=list(id$id, id$id)))
    if (any(duplicated(id$id))) stop("All id values must be unique")
    kmat <- diag(n+1) /2
    kmat[n+1,n+1]    <- 0 

    pdepth <- kindepth(id)
    mrow <- ifelse(id$mindex ==0, n+1, id$mindex)
    drow <- ifelse(id$findex ==0, n+1, id$findex)

    # Are there any MZ twins to worry about?
    if (!is.null(id$relation) && any(id$relation$code=="MZ twin")) {
        havemz <- TRUE
        temp <- which(id$relation$code=="MZ twin")
        ## drop=FALSE added in case only one MZ twin set
        mzmat <- as.matrix(id$relation[,c("indx1", "indx2")])[temp,,drop=FALSE]

        # any triples, quads, etc?
        if (any(table(mzmat) > 1)) { #yes there are
            # each group id will be min(member id)
            mzgrp <- 1:max(mzmat)  #each person a group
            indx <- sort(unique(as.vector(mzmat)))
            # The loop below will take k-1 iterations for a set labeled as
            #   1:2, 2:3, ...(k-1):k;  this is the worst case.
            while(1) {
                z1 <- mzgrp[mzmat[,1]]
                z2 <- mzgrp[mzmat[,2]]
                if (all(z1 == z2)) break
                mzgrp[indx] <- tapply(c(z1, z1, z2, z2), c(mzmat,mzmat), min)
            }
            # Now mzgrp = min person id for each person in a set
            matlist <- tapply(mzmat, mzgrp[mzmat], function(x) {
                x <- sort(unique(x))
                temp <- cbind(rep(x, each=length(x)), rep(x, length(x)))
                temp[temp[,1] != temp[,2],]
                })
            }
        else {  #no triples, easier case
            matlist <- tapply(mzmat, row(mzmat), function(x) 
                            matrix(x[c(1,2,2,1)],2), simplify=FALSE)
            }
        }
    else havemz <- FALSE

    for (depth in 1:max(pdepth)) {
        indx <- (1:n)[pdepth==depth]
        for (i in indx) {
            mom <- mrow[i]
            dad <- drow[i]
            kmat[i,]  <- kmat[,i] <- (kmat[mom,] + kmat[dad,])/2
            kmat[i,i] <- (1+ kmat[mom,dad])/2
            }
        if (havemz) {
            for (i in 1:length(matlist)) {
                temp <- matlist[[i]]
                kmat[temp] <- kmat[temp[1], temp[1]]
            }
        }
    }
    
    kmat <- kmat[1:n,1:n]
    dimnames(kmat) <- list(id$id, id$id)
    kmat
}    

oldkinship.pedigreeList <- function(id, ...) {
    famlist <- unique(id$famid)
    nfam <- length(famlist)
    matlist <- vector("list", nfam)
    idlist  <- vector("list", nfam) #the possibly reorderd list of id values
   
    for (i in 1:length(famlist)) {
        tped <- id[i]  #pedigree for this family
        temp <- try(oldkinship(tped, ...), silent=TRUE)
        if (class(temp)=="try-error") 
            stop(paste("In family", famlist[i], ":", temp))
        else matlist[[i]] <- as(forceSymmetric(temp), "dsCMatrix")
        idlist[[i]] <- tped$id
    }

    result <- bdiag(matlist)
    if (any(duplicated(id$id))) 
        temp <-paste(rep(famlist, sapply(idlist, length)),
                     unlist(idlist), sep='/') 
    else temp <- unlist(idlist)
        
    dimnames(result) <- list(temp, temp)
    result
}
