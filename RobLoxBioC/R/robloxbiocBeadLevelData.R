setMethod("robloxbioc", signature(x = "beadLevelData"),
    function(x, channelList = list(greenChannel), probeIDs = NULL, useSampleFac = FALSE, 
             sampleFac = NULL, weightNames = "wts", removeUnMappedProbes = TRUE,
             eps = NULL, eps.lower = 0, eps.upper = 0.05, steps = 3L, fsCor = TRUE, mad0 = 1e-4){
        BLData <- x
        arraynms <- sectionNames(BLData)
        output <- vector("list", length(channelList))
        if (useSampleFac) {
            if (is.null(sampleFac)) {
                if (!"SampleGroup" %in% names(BLData@sectionData)) {
                    cat("Could not determine sample factor from beadLevelData. Summarizing each section separately\n")
                    sList <- arraynms
                    sampleFac <- arraynms
                    newNames <- sList
                }else{
                    sampleFac <- BLData@sectionData$SampleGroup[, 1]
                    sList <- unique(sampleFac)
                    dupList <- which(duplicated(sampleFac))
                    if (any(dupList)) {
                        newNames <- strtrim(arraynms[-dupList], 12)
                    }else{ 
                        newNames <- strtrim(arraynms, 12)
                    }
                }
            }else{
                if (length(sampleFac) != length(arraynms)) {
                    cat("Length of specified sample factor did not match number of sections\n")
                    cat("length of sample factor: ", length(sampleFac), sampleFac, "\n")
                    cat("number of sections: ", length(arraynms), arraynms, "\n")
                    stop("Aborting summarization\n")
                }
                else {
                    sList <- unique(sampleFac)
                    dupList <- which(duplicated(sampleFac))
                    if (any(dupList)) {
                        newNames <- strtrim(arraynms[-dupList], 12)
                    }else{ 
                        newNames <- strtrim(arraynms, 12)
                    }
                }
            }
        }
        else {
            cat("No sample factor specified. Summarizing each section separately\n")
            sList <- arraynms
            sampleFac <- arraynms
            newNames <- arraynms
        }
        if (is.null(probeIDs)) {
            cat("Finding list of unique probes in beadLevelData\n")
            probeIDs <- beadarray:::uniqueProbeList(BLData)
            cat(length(probeIDs), " unique probeIDs found\n")
        }
        if (removeUnMappedProbes) {
            annoName <- annotation(BLData)
            if (!is.null(annoName)) {
                annoLoaded <- require(paste("illumina", annoName, ".db", sep = ""), character.only = TRUE)
                if (annoLoaded) {
                    mapEnv <- as.name(paste("illumina", annoName, "ARRAYADDRESS", sep = ""))
                    allMapped <- mappedkeys(revmap(eval(mapEnv)))
                    isMapped <- which(probeIDs %in% allMapped)
                    cat("Number of unmapped probes removed: ", length(probeIDs) - length(isMapped), "\n")
                    probeIDs <- probeIDs[isMapped]
                }
            }else{
                cat("Could not determine annotation for this beadLevelData object.\n")
            }
        }
        cNames <- unlist(lapply(channelList, function(x) x@name))
        if (any(duplicated(cNames))) {
            uNames <- unique(cNames)
            for (i in 1:length(uNames)) {
                sPos <- grep(uNames[i], cNames)
                if (length(sPos) > 1) {
                    for (j in 1:length(sPos)) {
                      cNames[sPos[j]] <- paste(cNames[sPos[j]], j, sep = ".")
                      warning("Duplicated channel names were found. Renaming...\n")
                    }
                }
            }
        }
        for (cNum in 1:length(channelList)) {
            template <- matrix(nrow = length(probeIDs), ncol = length(sList))
            if (length(channelList) == 1) { 
                newCols <- newNames
            }else{ 
                newCols <- paste(cNames[cNum], newNames, sep = ":")
            }
            output[[cNum]][["eMat"]] <- template
            colnames(output[[cNum]][["eMat"]]) <- newCols
            rownames(output[[cNum]][["eMat"]]) <- probeIDs
            output[[cNum]][["varMat"]] <- template
            colnames(output[[cNum]][["varMat"]]) <- newCols
            rownames(output[[cNum]][["varMat"]]) <- probeIDs
            output[[cNum]][["nObs"]] <- template
            colnames(output[[cNum]][["nObs"]]) <- newCols
            rownames(output[[cNum]][["nObs"]]) <- probeIDs
        }
        for (s in 1:length(sList)) {
            an <- which(sampleFac == sList[s])
            pIDs <- wts <- NULL
            values <- vector("list", length(channelList))
            for (i in an) {
                tmp <- BLData[[i]]
                pidCol <- grep("ProbeID", colnames(tmp))
                retainedBeads <- which(tmp[, pidCol] %in% probeIDs)
                tmp <- tmp[retainedBeads, ]
                wCol <- grep(weightNames, colnames(tmp))
                pIDs <- c(pIDs, tmp[, pidCol])
                if (length(wCol) == 0){ 
                    wts <- rep(1, length(pIDs))
                }else{ 
                    wts <- c(wts, tmp[, wCol])
                }
                for (ch in 1:length(channelList)) {
                    chName <- channelList[[ch]]@name
                    transFun <- channelList[[ch]]@transFun[[1]]
                    cat("Summarizing ", chName, " channel\n")
                    cat("Processing Array", i, "\n")
                    newVals <- transFun(BLData, array = i)[retainedBeads]
                    if (length(newVals) != nrow(tmp)) 
                      stop("Transformation function did not return correct number of values")
                    values[[ch]] <- c(values[[ch]], newVals)
                }
            }
            for (ch in 1:length(channelList)) {
                exprFun <- channelList[[ch]]@exprFun[[1]] #!
                varFun <- channelList[[ch]]@varFun[[1]]   #!
                values2 <- values[[ch]]
                naVals <- which(is.na(values2) | is.infinite(values2))
                pIDs2 <- pIDs
                wts2 <- wts
                if (length(naVals) > 0) {
                    values2 <- values2[-naVals]
                    pIDs2 <- pIDs[-naVals]
                    wts2 <- wts[-naVals]
                }
                pOrder <- order(pIDs2)
                pIDs2 <- pIDs2[pOrder]
                values2 <- values2[pOrder]
                wts2 <- wts2[pOrder]
                if (any(wts2 == 0)) {
                    values2 <- values2[-which(wts2 == 0)]
                    pIDs2 <- pIDs2[-which(wts2 == 0)]
                    wts2 <- wts2[-which(wts2 == 0)]
                }
                ## spezieller Teil für rmx
                tmp <- split(wts2 * values2, pIDs2)
                pMap <- match(names(tmp), probeIDs)
                cat("Using rmx\n")
                res.rmx <- rmxBeadSummary(x = tmp, eps = eps, eps.lower = eps.lower, eps.upper= eps.upper, 
                                          steps = steps, fsCor = fsCor, mad0 = mad0)
                output[[ch]][["eMat"]][pMap, s] <- res.rmx$eMat
                output[[ch]][["varMat"]][pMap, s] <- res.rmx$varMat
                output[[ch]][["nObs"]][pMap, s] <- res.rmx$nObs
            }
        }
        cat("Making summary object\n")
        eMat <- output[[1]][["eMat"]]
        varMat <- output[[1]][["varMat"]]
        nObs <- output[[1]][["nObs"]]
        if (length(output) > 1) {
            for (i in 2:length(output)) {
                eMat <- cbind(eMat, output[[i]][["eMat"]])
                varMat <- cbind(varMat, output[[i]][["varMat"]])
                nObs <- cbind(nObs, output[[i]][["nObs"]])
            }
        }
        channelFac <- NULL
        for (i in 1:length(channelList)) {
            newfac <- cNames[i]
            channelFac <- c(channelFac, rep(newfac, length(sList)))
        }
        BSData <- new("ExpressionSetIllumina")
        annoName <- annotation(BLData)
        if (!is.null(annoName)) {
            annoLoaded <- require(paste("illumina", annoName, ".db", sep = ""), character.only = TRUE)
            if (annoLoaded) {
                mapEnv <- as.name(paste("illumina", annoName, "ARRAYADDRESS", sep = ""))
                IlluminaIDs <- as.character(unlist(AnnotationDbi::mget(as.character(probeIDs), revmap(eval(mapEnv)), ifnotfound = NA)))
                rownames(eMat) <- rownames(varMat) <- rownames(nObs) <- as.character(IlluminaIDs)
                status <- rep("Unknown", length(probeIDs))
                annoPkg <- paste("illumina", annoName, ".db", sep = "")
                annoVers <- packageDescription(annoPkg, field = "Version")
                message(paste("Annotating control probes using package ", annoPkg, " Version:", annoVers, "\n", sep = ""))
                mapEnv <- as.name(paste("illumina", annoName, "REPORTERGROUPNAME", sep = ""))
                t <- try(eval(mapEnv), silent = TRUE)
                if (class(t) == "try-error") {
                    message(paste("Could not find a REPORTERGROUPNAME mapping in annotation package ", annoPkg, 
                                  ". Perhaps it needs updating?", sep = ""))
                }
                else {
                    status[which(!is.na(IlluminaIDs))] <- unlist(AnnotationDbi::mget(IlluminaIDs[which(!is.na(IlluminaIDs))], 
                                                                      eval(mapEnv), ifnotfound = NA))
                    status[which(is.na(status))] <- "regular"
                }
                featureData(BSData) <- new("AnnotatedDataFrame", data = data.frame(ArrayAddressID = probeIDs, 
                    IlluminaID = IlluminaIDs, Status = status, row.names = IlluminaIDs))
                BSData@annotation <- annoName
            }
        }
        else {
            cat("Could not map ArrayAddressIDs: No annotation specified\n")
            featureData(BSData) <- new("AnnotatedDataFrame", data = data.frame(ProbeID = probeIDs, 
                                                                               row.names = probeIDs))
        }
        assayData(BSData) <- assayDataNew(exprs = eMat, se.exprs = varMat, 
                                          nObservations = nObs, storage.mode = "list")
        sampInfo <- beadarray:::sampleSheet(BLData)
        if (!is.null(sampInfo)) {
            expIDs <- paste(sampInfo$Sentrix_ID, sampInfo$Sentrix_Position, sep = "_")
            sampInfo <- sampInfo[sapply(newNames, function(x) grep(strtrim(x, 12), expIDs)), ]
            p <- new("AnnotatedDataFrame", data.frame(sampInfo, row.names = newNames))
        }
        else p <- new("AnnotatedDataFrame", data.frame(sampleID = newNames, SampleFac = unique(sampleFac), 
                                                       row.names = newNames))
        phenoData(BSData) <- p
        qcNames <- names(BLData@sectionData)
        qcNames <- setdiff(qcNames, "Targets")
        if (length(qcNames) > 0) {
            QC <- BLData@sectionData[[qcNames[[1]]]]
            if (length(qcNames > 1)) {
                for (i in 2:length(qcNames)) {
                    QC <- cbind(QC, BLData@sectionData[[qcNames[i]]])
                }
            }
            QC <- new("AnnotatedDataFrame", data.frame(QC, row.names = arraynms))
            BSData@QC <- QC
        }
        BSData@channelData <- list(channelFac, channelList)
        BSData
    })
## computation of bead summaries via robloxbioc for "matrix"
rmxBeadSummary <- function(x, eps, eps.lower, eps.upper, steps, fsCor, mad0){
    nObs <- sapply(x, length)
    noBeads <- as.integer(names(table(nObs)))
    varMat <- eMat <- numeric(length(nObs))
    for(i in seq(along = noBeads)){
        index <- nObs == noBeads[i]
        if(noBeads[i] == 1){
          eMat[index] <- as.vector(unlist(x[index]))
          varMat[index] <- mad0
        }else{
            temp <- t(sapply(x[index], rbind))
            rmx <- robloxbioc(temp, eps = eps, eps.lower = eps.lower, eps.upper = eps.upper, 
                              steps = steps, fsCor = fsCor, mad0 = mad0)
            eMat[index] <- rmx[,"mean"]
            varMat[index] <- rmx[,"sd"]
        }
    }
    list(eMat = eMat, varMat = varMat, nObs = nObs)
}
