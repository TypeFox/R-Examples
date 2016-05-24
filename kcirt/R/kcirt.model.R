kcirt.model <-
function(constructMap.ls, qTypes, data=NULL, Y=NULL, mu=0, mxLambda=NULL, covEta=1, covShocks=1, deltaType=1) {
    
    nBlocks <- length(constructMap.ls) ; nBlocks
    
    constructs <- unique( unlist(constructMap.ls) ) ; constructs
    
    nuc <- length(constructs) ; nuc
    
    ns <- unlist( lapply(constructMap.ls, length) )
    nks <- as.integer((ns-1)*ns/2)
    
    
    ################################################ create mx Delta
    if(deltaType==1) {
        mxDelta <- matrix(0, sum(nks), sum(ns))
        ndx.startRow <- 0
        ndx.startColumn <- 0
        for(ib in 1:nBlocks) {
            ii.vec <- rep(NA, nks[ib])
            jj.vec <- rep(NA, nks[ib])
            kk <- 0
            for(ii in 2:ns[ib]) {
                for(jj in 1:(ii-1)) {
                    kk <- kk+1
                    ii.vec[kk] <- ii
                    jj.vec[kk] <- jj
                }
            }
            for(kk in 1:nks[ib]) {
                mxDelta[ ndx.startRow+kk, ndx.startColumn + ii.vec[kk] ] <- 1
                mxDelta[ ndx.startRow+kk, ndx.startColumn + jj.vec[kk] ] <- -1
            }
            ndx.startRow <- ndx.startRow + nks[ib]
            ndx.startColumn <- ndx.startColumn + ns[ib]
        }
        
        
        #ii.vec <- c(2, 3, 3, 4, 4, 4)
        #jj.vec <- c(1, 1, 2, 1, 2, 3)
    }
    
    ########################### Like OLIVERAS BROWN
    if(deltaType==2) {
        mxDelta <- matrix(0, sum(nks), sum(ns))
        ndx.startRow <- 0
        ndx.startColumn <- 0
        for (ib in 1:nBlocks) {
            ii.vec <- rep(NA, nks[ib])
            jj.vec <- rep(NA, nks[ib])
            kk <- 0
            for (ii in 1:(ns[ib]-1)) {
                for (jj in (ii+1):ns[ib]) {
                    kk <- kk + 1
                    ii.vec[kk] <- ii
                    jj.vec[kk] <- jj
                }
            }
            for (kk in 1:nks[ib]) {
                mxDelta[ndx.startRow + kk, ndx.startColumn + ii.vec[kk]] <- 1
                mxDelta[ndx.startRow + kk, ndx.startColumn + jj.vec[kk]] <- -1
            }
            ndx.startRow <- ndx.startRow + nks[ib]
            ndx.startColumn <- ndx.startColumn + ns[ib]
        }
         
    }
    
    ################################################ create mx Slot
    
    constructKey <- sort(unique(unlist(constructMap.ls))) ; constructKey
    
    mxSlot <- matrix(0, sum(ns), nuc)
    #xzero.vec <- rep(0, nuc)
    ndx.startRow <- 0
    for(ib in 1:nBlocks) {
        thisCmap <- constructMap.ls[[ib]] ; thisCmap
        thisCmap <- match(thisCmap, constructKey) ; thisCmap
        for(i in 1:length(thisCmap)) {
            ndx.startRow <- ndx.startRow + 1
            mxSlot[ ndx.startRow, thisCmap[i] ] <- 1
        }
    }
    
    
    
    
    ################################################ create generic mx Lambda
    
    if(is.null(mxLambda)) {
        mxLambda <- diag( rep_len( c(1,-1), sum(ns) ) )
    } else {
        if(is.vector(mxLambda)) {
            mxLambda <- diag( mxLambda, sum(ns) )
        }
    }
    
    
    if( !is.matrix(covEta) ) {
        covEta <- diag(covEta, nuc)
    }
    
    if( !is.matrix(covShocks) ) {
        covShocks <- diag(covShocks, sum(ns))
    }
    
    
    mu <- rep_len(mu, sum(ns))
    
    
    ############################################### create crosstalk Lambda Info
    Qid <- list()
    for(ii in 1:nBlocks) {
        Qid[[ii]] <- rep(ii, ns[ii])
    }
    
    
    mxLambdaCTinfo <- matrix("S", sum(ns), sum(ns))
    xcmapVec <- unlist(constructMap.ls) ; xcmapVec
    xqidVec <- unlist(Qid)
    
    for(ii in 2:sum(ns)) {
        for(jj in 1:(ii-1)) {
            
            if(xcmapVec[ii] == xcmapVec[jj]) {
                sameConstruct <- "T"
            } else {
                sameConstruct <- "F"
            }
            
            if(xqidVec[ii] == xqidVec[jj]) {
                sameBlock <- "W"
            } else {
                sameBlock <- "B"
            }
            keyStr <- paste(sameBlock, sameConstruct, sep="")
            
            mxLambdaCTinfo[ii, jj] <- keyStr
            mxLambdaCTinfo[jj, ii] <- keyStr
        }
    }
    
    covStochastic <-  mxDelta %*% covShocks %*% t(mxDelta)
    mxDLS <- mxDelta %*% mxLambda %*% mxSlot
    mxSigma <- mxDLS %*% covEta %*% t(mxDLS)   +  covStochastic
    
    
    
    if(!is.null(data)) {
        Y.created <- ikcirt.data2Y(mxData=data, mxDelta=mxDelta)
    }
    if(!is.null(Y)) {
        data.created <- ikcirt.Y2data(Y=Y, mxDelta=mxDelta, ns=ns)
    }
    
    if(exists("Y.created") & !is.null(Y)) {
        ### CHECK THAT THEY ARE THE SAME
    }
    
    if(exists("data.created") & !is.null(data)) {
        ### CHECK THAT THEY ARE THE SAME
    }
    
    if(exists("Y.created")) {
        Y.out <- Y.created
    } else {
        Y.out <- Y
    }
    if(exists("data.created")) {
        data.out <- data.created
    } else {
        data.out <- data
    }
    
    if(is.matrix(Y.out)) {
        Z <- 2*Y.out-1
        Yisna <- is.na(Y.out)
    } else {
        Z <- NULL
        Yisna <- NULL
    }
    
    out.ls <- list(
    "constructMap.ls" = constructMap.ls,
    "qTypes" = qTypes,
    "mxData" = data.out,
    "Y" = Y.out,
    "Z"=Z,
    "Yisna"=Yisna,
    "nBlocks" = nBlocks,
    "nuc" = nuc,
    "ns" = ns,
    "nks" = nks,
    "mxDelta" = mxDelta,
    "mxSlot" = mxSlot,
    "mu" = mu,
    "mxLambda" = mxLambda,
    "covEta"=covEta,
    "covShocks"=covShocks,
    "covStochastic"=covStochastic,
    "mxSigma"=mxSigma,
    "Qid" = Qid,
    "mxLambdaCTinfo" = mxLambdaCTinfo,
    "constructSlotKey" = constructKey
    )
    
    class(out.ls) <- "kcube.irt.model"
    
    return(out.ls)
}
