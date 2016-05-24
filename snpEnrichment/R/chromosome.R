######################################################################
########################## Class Chromosome ##########################
############################## Creation ##############################
######################################################################


setClass(
    Class = "Chromosome",
    representation = representation(
        Data = "data.frame",
        LD = "character",
        eSNP = "EnrichSNP",
        xSNP = "EnrichSNP"
    ),
    prototype = prototype(
        Data = data.frame(),
        LD = character(),
        eSNP = enrichSNP(),
        xSNP = enrichSNP()
    )
)


setMethod(f = "chromosome", signature = "ANY", definition = function (Data, LD, eSNP, xSNP) {
    if (missing(Data)) {
        Data <- data.frame()
        if (missing(eSNP)) {
            eSNP <- enrichSNP()
            xSNP <- enrichSNP()
        } else {}
    } else {
        if (missing(eSNP)) {
            eSNP <- enrichSNP(List = Data[Data[, "eSNP"] == 1, "SNP"])
            xSNP <- enrichSNP(List = Data[Data[, "xSNP"] == 1, "SNP"])
        } else {}
    }
    if (missing(LD)) {
        LD <- character()
    } else {}
    return(new("Chromosome", Data = Data, LD = LD, eSNP = eSNP, xSNP = xSNP))
})


setMethod(f = "is.chromosome", signature = "ANY", definition = function (object) {
    if (length(object)>1) {
        return(sapply(object, is.chromosome))
    } else {
        if (class(object) == "Chromosome") {
            return(TRUE)
        } else {
            return(FALSE)
        }
    }
})


setMethod(f = "print", signature = "Chromosome", definition = function (x, type = c("eSNP", "xSNP")) {
    if (missing(x)) {
        stop('[Chromosome:print] "x" is missing.', call. = FALSE)
    } else {}
    if (is.null(type) | any(!type%in%c("eSNP", "xSNP"))) {
        stop('[Chromosome:print] "type" must be: "eSNP" and/or "xSNP".', call. = FALSE)
    } else {}

    res <- list()
    for (iType in type) {
        resTmp <- print(x[iType])
        colnames(resTmp) <- c("EnrichmentRatio", "Z", "PValue", "nSample", "TotalSNP", iType)
        rownames(resTmp) <-  paste("Chrom", iType, sep = ":")
        res[[iType]] <- resTmp
    }
    if (length(type)==1) {
        return(res[[1]])
    } else {
        res <- do.call("rbind", res)
        rownames(res) <- paste("Chrom", type, sep = ":")
        return(res)
    }
    return(res)
})


.Chromosome.show <- function (object) {
    cat("  ~ Data :", paste0("(", paste(dim(object@Data), collapse = "x"), ")"))
        nrowShow <- seq(min(5 , nrow(object@Data)))
        ncolShow <- seq(min(10 , ncol(object@Data)))
        if (nrow(object@Data) == 0) {
            cat(" NA")
        } else {
            cat("\n")
            resFormat <- cbind(c("", rownames(object@Data[nrowShow, ncolShow])), rbind(colnames(object@Data[nrowShow, ncolShow]), apply(object@Data[nrowShow, ncolShow], 2, as.character)))
            resFormat <- rbind(resFormat, ".....")
            cat(paste("     ", apply(apply(matrix(paste0(" ", resFormat, " "), nrow = nrow(resFormat)), 2, format, justify = "centre"), 1, paste, collapse = ""), "\n", sep = "", collapse = ""))
        }
    cat("\n  ~ LD :")
        if (length(object@LD) == 0) {
            cat(" NA")
        } else {
            cat("\n")
            if (length(object@LD) > 5) {
                tmpLD <- object@LD[seq(5)]
            } else {
                tmpLD <- object@LD
            }
            resFormat <- cbind(c("SNP1", "SNP2"), rbind(names(tmpLD), tmpLD))
            resFormat <- cbind(resFormat, "...")
            cat(paste("     ", apply(apply(matrix(paste0(" ", resFormat, " "), nrow = nrow(resFormat)), 2, format, justify = "centre"), 1, paste, collapse = ""), "\n", sep = "", collapse = ""))
        }
    cat("\n  ~ eSNP :\n")
        .EnrichSNP.show(object@eSNP)
    cat("\n  ~ xSNP :\n")
        .EnrichSNP.show(object@xSNP)
    cat("\n")
    return(invisible(object))
}
setMethod(f = "show", signature = "Chromosome", definition = function (object) {cat("     ~~~ Class:", class(object), "~~~\n"); .Chromosome.show(object); return(invisible(object))})


setMethod(f = "[", signature = "Chromosome", definition = function (x, i, j, drop) {
    switch(EXPR = i,
        "Data" = {return(x@Data)},
        "LD" = {return(x@LD)},
        "eSNP" = {return(x@eSNP)},
        "xSNP" = {return(x@xSNP)},
        stop('[Chromosome:get] ', i, ' is not a "Chromosome" slot.', call. = FALSE)
    )
})


setMethod(f = "[<-", signature = "Chromosome", definition = function (x, i, j, value) {
    switch(EXPR = i,
        "Data" = {x@Data <- value},
        "LD" = {x@LD <- value},
        "eSNP" = {x@eSNP <- value},
        "xSNP" = {x@xSNP <- value},
        stop('[Chromosome:get] ', i, ' is not a "Chromosome" slot.', call. = FALSE)
    )
    validObject(x)
    return(invisible(x))
})


setMethod(f = "computeER", signature = "Chromosome", definition = function (object, sigThresh = 0.05, mc.cores = 1) {
    if (!missing(object)) {
        object@Chromosomes <- mclapply2(object@Chromosomes, mc.cores = mc.cores, function (chr) {
            data <- chr@Data
            chrLD <- length(chr@LD)
            for (type in c("eSNP", "xSNP")) {
                if (!(chrLD == 0 & type == "xSNP")) {
                    snpEnrich <- table(factor(chr@Data[, "PVALUE"]<sigThresh, levels = c(FALSE, TRUE)), factor(chr@Data[, type], levels = c(0, 1)))
                    colnames(snpEnrich) <- c("otherSNP", type)
                    rownames(snpEnrich) <- eval(parse(text = paste0('c("P>=', sigThresh, '", "P<', sigThresh, '")')))
                    chr[type] <- enrichSNP(List = chr[type]@List, Table = unclass(snpEnrich), EnrichmentRatio = .enrichmentRatio(snpEnrich))
                } else {}
            }
            return(chr)
        })
        return(object)
    } else {
        stop('[Chromosome:computeER] "Chromosome" object is required.', call. = FALSE)
    }
})


setMethod(f = "doLDblock", signature = "Chromosome", definition = function (object, mc.cores = 1) {
    if (!missing(object)) {
        data <- object@Data

        nbCORES <- mc.cores
        nbCores <- max(1, round((nbCORES-22)/22))
        chrLD <- object@LD

        byBlock <- split(chrLD, names(chrLD))
        byBlock <- unique(mclapply2(byBlock, mc.cores = nbCores, function (i) {
            names(i) <- NULL
            return(i)
        }))
        LDblockTmp <- mclapply2(seq_along(byBlock), mc.cores = nbCores, function (jBlock) {
            isIn <- which(data[, "SNP"] %in% byBlock[[jBlock]])
            if (length(isIn) > 0) {
                return(range(data[which(data[, "SNP"] %in% byBlock[[jBlock]]), "POS"]))
            } else {
                return(NA)
            }
        })
        LDblock <- do.call("rbind", unique(LDblockTmp))
        LDblock <- na.exclude(LDblock)
        rm(LDblockTmp)

        names(byBlock) <- NULL
        LDblock <- LDblock[order(LDblock[, 1]), ]
        rm(chrLD, byBlock)

        blockLim <- NULL
        for (iBlock in seq(nrow(LDblock))) {
            if (iBlock == 1) {
                POS <- LDblock[iBlock, ]
                blockLim <- rbind(blockLim, POS)
                jBlock <- 1
            } else {
                POS <- LDblock[iBlock, ]
                if (max(blockLim[nrow(blockLim), ])<min(POS)) {
                    blockLim <- rbind(blockLim, POS)
                } else {
                    blockLim[nrow(blockLim), ] <- range(blockLim[nrow(blockLim), ], POS)
                }
                if (iBlock == nrow(LDblock)) {
                    iBlock <- iBlock + 1
                    if (max(blockLim[nrow(blockLim), ])<min(POS)) {
                        blockLim <- rbind(blockLim, POS)
                    } else {
                        blockLim[nrow(blockLim), ] <- range(blockLim[nrow(blockLim), ], POS)
                    }
                } else {}
            }
        }
        rm(LDblock)

        blockLim <- cbind(blockLim, seq(nrow(blockLim)))
        colnames(blockLim) <- c("MIN", "MAX", "IDBLOCK")
        rownames(blockLim) <- seq(nrow(blockLim))
        blockLim <- cbind(blockLim, LENGTH = NA)
        resParallel <- mclapply2(seq(nrow(blockLim)), mc.cores = nbCores, function (li) {
            blockLim[li, "LENGTH"] <- as.integer(blockLim[li, "MAX"])-as.integer(blockLim[li, "MIN"])
            return(blockLim[li, ])
        })
        blockLim <- do.call("rbind", resParallel)
        rm(resParallel)

        data <- data[order(data[, "POS"]), ]
        data[, c("MIN", "MAX", "IDBLOCK", "LENGTH", "MAFmedian")] <- as.numeric(NA)
        tmpChr <- mclapply2(seq(nrow(blockLim)), mc.cores = nbCores, function (i) {
            m <- blockLim[i, ]
            interv <- seq(from = which(data[, "POS"] == m["MIN"]), to = which(data[, "POS"] == m["MAX"]))
            interv
            data[interv, c("MIN", "MAX", "IDBLOCK", "LENGTH")] <- matrix(rep(m, length(interv)), nrow = length(interv), byrow = TRUE)
            data[which(data[, "IDBLOCK"]%in%m["IDBLOCK"]), "MAFmedian"] <- median(data[data[, "IDBLOCK"]%in%m["IDBLOCK"], "MAF"])
            return(data[interv, ])
        })

        dataTmp <- do.call("rbind", tmpChr)
        rm(tmpChr, blockLim)

        missingData <- data[!data[, "SNP"]%in%dataTmp[, "SNP"], ]
        maxIDBLOCK <- max(dataTmp[, "IDBLOCK"])
        for (iRow in seq(nrow(missingData))) {
            missingData[iRow, "MIN"] <- missingData[iRow, "POS"]
            missingData[iRow, "MAX"] <- missingData[iRow, "POS"]
            missingData[iRow, "LENGTH"] <- missingData[iRow, "MAX"] - missingData[iRow, "MIN"]
            missingData[iRow, "MAFmedian"] <- missingData[iRow, "MAF"]
            missingData[iRow, "IDBLOCK"] <- maxIDBLOCK + iRow
        }

        data <- rbind(missingData, dataTmp)
        data <- data[!is.na(data[, "SNP"]), ]
        data <- data[order(data[, "POS"]), ]
        rownames(data) <- data[, "SNP"]
        object@Data <- data
        return(object)
    } else {
        stop('[Chromosome:doLDblock] "Chromosome" Object is required.', call. = FALSE)
    }
})


setMethod(f = "reSample", signature = "Chromosome", definition = function (object, nSample = 100, empiricPvalue = TRUE, sigThresh = 0.05, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1) {
    if (!missing(object)) {
        if (nSample<10) {
            nSample = 10
            warning('[Enrichment:reSample] nSample was increased to 10.', call. = FALSE)
        } else {}
        result <- .reSample(object = object, nSample = nSample, empiricPvalue = empiricPvalue, sigThresh = sigThresh, MAFpool = MAFpool, mc.cores = mc.cores)
        return(result)
    } else {
        stop('[Enrichment:reSample] "Enrichment" object is required.', call. = FALSE)
    }
})


setMethod(f = "reset", signature = "Chromosome", definition = function (object, i) {
    switch(EXPR = i,
        "Data" = {
            object@Data <- data.frame()
        },
        "LD" = {
            object@LD <- character()
        },
        "eSNP" = {
            object@eSNP <- enrichSNP()
        },
        "xSNP" = {
            object@xSNP <- enrichSNP()
        },
        "List" = {
            for (type in c("eSNP", "xSNP")) {
                object[type]@List <- reset(object[type], "List")
            }
        },
        "Table" = {
            for (type in c("eSNP", "xSNP")) {
                object[type]@Table <- matrix(0, ncol = 2, nrow = 2)
            }
        },
        "EnrichmentRatio" = {
            for (type in c("eSNP", "xSNP")) {
                object[type]@EnrichmentRatio <- numeric()
            }
        },
        "Z" = {
            for (type in c("eSNP", "xSNP")) {
                object[type]@Z <- numeric()
            }
        },
        "PValue" = {
            for (type in c("eSNP", "xSNP")) {
                object[type]@PValue <- numeric()
            }
        },
        "Resampling" = {
            for (type in c("eSNP", "xSNP")) {
                object[type]@Resampling <- matrix(0, ncol = 5, nrow = 0)
            }
        },
        stop('[Enrichment:reset] ', i, ' is not a "Enrichment" slot.', call. = FALSE)
    )
    return(object)
})
