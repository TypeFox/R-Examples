fcrosChrSegment = function(chrData, nd = 10) {
    ## get fcros segmentation values
    cdata <- chrData$f.call
    idx_d <- which(cdata == -1)
    idx_g <- which(cdata == 1)
    idStart <- sort(c(idx_d, idx_g))
    idEnd <- idStart
    segProba <- chrData$f.value[idStart]
    positions <- chrData$Position
    lBound <- positions[idStart]
    uBound <- lBound
    nbSeg <- length(lBound)
    L2R <- chrData$f.L2R
    segVal <- L2R[idStart]
    ndata <- length(L2R)
    sigma <- mad((L2R[-1] - L2R[-ndata])/sqrt(2))
    dm <- (positions[ndata] - positions[1]) / (ndata-1)

    ############################## first merge of segments
    j <- 1
    while(j < nbSeg) {
        j1 <- idEnd[j]
        j2 <- idStart[j+1]
        nb <- j2-j1+1
        seg_d <- 0
        dnb <- nb*dm     # distance between the 2 segments
        d1 <- dnb
        if (nb < nd) d1 <- lBound[j+1] - uBound[j]
        if (d1 < dnb)  seg_d <- 1
        seg_v <- 0
        v1 <- ((cdata[j1] > 0) && (cdata[j2] > 0)) # amplification
        v2 <- ((cdata[j1] < 0) && (cdata[j2] < 0)) # deletion
        if (v1 || v2) seg_v <- 1

        if (seg_d && seg_v) {	          # now merge the j and j+1 segments
           uBound[j] <- uBound[j+1]	  # modify upper bounds
           idEnd[j] <- idEnd[j+1]
           ns <- idEnd[j]-idStart[j]+1
           xbar <- mean(L2R[idStart[j]:idEnd[j]])
           segVal[j] <- xbar              # modify mean val
           segStat <- sqrt(ns)*xbar/sigma # modify stat val
           segProba[j] <- pnorm(segStat, 0, 1)   # central limit theorem
           for (k in (j+1):(nbSeg-1)) {   # modify the other segments infos
               lBound[k] <- lBound[k+1]
               uBound[k] <- uBound[k+1]
               segVal[k] <- segVal[k+1]
               segProba[k] <- segProba[k+1]
               idStart[k] <- idStart[k+1]
               idEnd[k] <- idEnd[k+1]
           }
           nbSeg <- nbSeg - 1
        } else { j <- j+1  }
    }
    lBound <- lBound[1:nbSeg]
    uBound <- uBound[1:nbSeg]
    segVal <- segVal[1:nbSeg]
    segProba <- segProba[1:nbSeg]
    idStart <- idStart[1:nbSeg]
    idEnd <- idEnd[1:nbSeg]

    ############################## single probe segment treatment
    ## check if last probe is single region
    tmp <- idEnd[nbSeg] - idStart[nbSeg]
    dtmp <- positions[idStart[nbSeg]] - positions[idStart[nbSeg]-1]
    if ((!tmp) && (dtmp < dm)) nbSeg <- nbSeg-1

    ## deletion of region with single probe not in sparse region
    j <- 1
    while (j < nbSeg) {
         tmp <- idEnd[j] - idStart[j]
         if (!tmp) {  # check if the single probe should be kept
            if (idStart[j] == 1) {
               ## the first probe is a region
               dtmp <- positions[idEnd[j]+1] - positions[idEnd[j]]
            } else {
               dbefor <- positions[idStart[j]] - positions[idStart[j]-1]
               dafter <- positions[idEnd[j]+1] - positions[idEnd[j]]
               dtmp <- min(dbefor, dafter)
            }
            if (dtmp < dm) {   # not in a sparse region, delete it
               cdata[idStart[j]] <- 0;
               for (k in j:nbSeg) {
                   # modify the other segments infos
                   lBound[k] <- lBound[k+1]
                   uBound[k] <- uBound[k+1]
                   segVal[k] <- segVal[k+1]
                   segProba[k] <- segProba[k+1]
                   idStart[k] <- idStart[k+1]
                   idEnd[k] <- idEnd[k+1]
               }
               nbSeg <- nbSeg-1
            }  else { j <- j+1} #  in a sparse region, keep it
         } else {j <- j+1 }
    }
    lBound <- lBound[1:nbSeg]
    uBound <- uBound[1:nbSeg]
    segVal <- segVal[1:nbSeg]
    segProba <- segProba[1:nbSeg]
    idStart <- idStart[1:nbSeg]
    idEnd <- idEnd[1:nbSeg]

    ############################## second merge of segments
    j <- 1
    while(j < nbSeg) {
        j1 <- idEnd[j]
        j2 <- idStart[j+1]
        nb <- j2-j1+1
        seg_d <- 0
        dnb <- nb*dm     # distance between the 2 segments
        d1 <- dnb
        if (nb < nd) d1 <- lBound[j+1] - uBound[j]
        if (d1 < dnb)  seg_d <- 1
        seg_v <- 0
        v1 <- ((cdata[j1] > 0) && (cdata[j2] > 0)) # amplification
        v2 <- ((cdata[j1] < 0) && (cdata[j2] < 0)) # deletion
        if (v1 || v2) seg_v <- 1

        if (seg_d && seg_v) {	          # now merge the j and j+1 segments
           uBound[j] <- uBound[j+1]	  # modify upper bounds
           idEnd[j] <- idEnd[j+1]
           ns <- idEnd[j+1]-idStart[j]+1
           xbar <- mean(L2R[idStart[j]:idEnd[j]])
           segVal[j] <- xbar              # modify mean val
           segStat <- sqrt(ns)*xbar/sigma # modify stat val
           segProba[j] <- pnorm(segStat, 0, 1)   # central limit theorem
           for (k in (j+1):(nbSeg-1)) {   # modify the other segments infos
               lBound[k] <- lBound[k+1]
               uBound[k] <- uBound[k+1]
               segVal[k] <- segVal[k+1]
               segProba[k] <- segProba[k+1]
               idStart[k] <- idStart[k+1]
               idEnd[k] <- idEnd[k+1]
           }
           nbSeg <- nbSeg - 1
        } else { j <- j+1  }
    }
    lBound <- lBound[1:nbSeg]
    uBound <- uBound[1:nbSeg]
    segVal <- segVal[1:nbSeg]
    segProba <- segProba[1:nbSeg]
    idStart <- idStart[1:nbSeg]
    idEnd <- idEnd[1:nbSeg]

    segData <- idStart[1:nbSeg]
    segData <- cbind(segData, idEnd[1:nbSeg])
    segData <- cbind(segData, lBound[1:nbSeg])
    segData <- cbind(segData, uBound[1:nbSeg])
    segData <- cbind(segData, segVal[1:nbSeg])
    segData <- cbind(segData, segProba[1:nbSeg])

    colnames(segData) <- c("idStart", "idEnd", "lBound", "uBound", "segL2R", "segProba")

    return(segData)
}
