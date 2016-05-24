setGeneric(name = "enrichSNP", def = function (List, Table, EnrichmentRatio, Z, PValue, Resampling) {standardGeneric("enrichSNP")})
setGeneric(name = "chromosome", def = function (Data, LD, eSNP, xSNP) {standardGeneric("chromosome")})
setGeneric(name = "enrichment", def = function (Loss, Call, eSNP, xSNP, Chromosomes) {standardGeneric("enrichment")})
setGeneric(name = "doLDblock", def = function (object, mc.cores = 1) {standardGeneric("doLDblock")})
setGeneric(name = "excludeSNP", def = function (object, excludeFile, mc.cores = 1) {standardGeneric("excludeSNP")})
setGeneric(name = "computeER", def = function (object, sigThresh = 0.05, mc.cores = 1) {standardGeneric("computeER")})
setGeneric(name = "reSample", def = function (object, nSample = 100, empiricPvalue = TRUE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = TRUE, ...) {standardGeneric("reSample")})
setGeneric(name = "compareEnrichment", def = function (object.x, object.y, pattern = "Chrom", nSample = 100, empiricPvalue = TRUE, mc.cores = 1, onlyGenome = TRUE) {standardGeneric("compareEnrichment")})
setGeneric(name = "is.enrichment", def = function (object) {standardGeneric("is.enrichment")})
setGeneric(name = "is.chromosome", def = function (object) {standardGeneric("is.chromosome")})
setGeneric(name = "is.EnrichSNP", def = function (object) {standardGeneric("is.EnrichSNP")})
setGeneric(name = "reset", def = function (object, i) {standardGeneric("reset")})
setGeneric(name = "getEnrichSNP", def = function (object, type = "eSNP") {standardGeneric("getEnrichSNP")})

setMethod(f = "reSample", signature = "ANY", definition = function (object, nSample, sigThresh, MAFpool, mc.cores) {
    if (!(is.enrichment(object) & is.chromosome(object))) {
        stop('[Method:reSample] not available for "', class(object), '" object.', call. = FALSE)
    } else {}
})
setMethod(f = "excludeSNP", signature = "ANY", definition = function (object, excludeFile, mc.cores = 1) {
    if (!is.enrichment(object)) {
        stop('[Method:excludeSNP] not available for "', class(object), '" object.', call. = FALSE)
    } else {}
})
setMethod(f = "reset", signature = "ANY", definition = function (object, i) {
    if (!(is.enrichment(object) & is.chromosome(object))) {
        stop('[Method:reset] not available for "', class(object), '" object.', call. = FALSE)
    } else {}
})
setMethod(f = "getEnrichSNP", signature = "ANY", definition = function (object, type = "eSNP") {
    if (!(is.enrichment(object))) {
        stop('[Method:getEnrichSNP] not available for "', class(object), '" object.', call. = FALSE)
    } else {}
})


.verbose <- function (expr) {return(invisible(capture.output(expr)))}


GC <- function (verbose = getOption("verbose"), reset = FALSE) {
    while (!identical(gc(verbose, reset)[, 4], gc(verbose, reset)[, 4])) {}
    return(gc(verbose, reset))
}


.checkFilePath <- function (path) {
    # END <- unlist(regmatches(path, regexec(".$", path)))
    # START <- unlist(regmatches(path, regexec("^.", path)))
    # if (START != "/") {
        # path <- paste0("/", path)
    # } else {}
    # if (END != "/") {
        # path <- paste0(path, "/")
    # } else {}
    path <- gsub("/*$", "/", path)
    return(path)
}


maxCores <- function (mc.cores = 1) {
    if (Sys.info()[["sysname"]] == "Linux") {
        nbCores <- detectCores()
        mc.cores.old <- mc.cores
        if (file.exists("/proc/meminfo")) {
            memInfo <- readLines("/proc/meminfo")
            sysMemFree <- memInfo[grep('^MemFree:', memInfo)]
            sysMemCached <- memInfo[grep('^Cached:', memInfo)]
            sysMemAvailable <- 0.95*(as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", sysMemFree)) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", sysMemCached)))
            sysProc <- as.numeric(unlist(strsplit(system(paste("ps v", Sys.getpid()), intern = TRUE)[2], " +"))[8])
            mc.cores <- max(min(as.numeric(mc.cores), floor(sysMemAvailable/sysProc)), 1)
            if (mc.cores > nbCores) {
                mc.cores <- nbCores
            } else {}
            if (mc.cores != mc.cores.old) {
                warning(paste0('To avoid memory overload "mc.cores" was decreased to "', mc.cores, '".'), call. = FALSE)
            } else {}
        } else {
            mc.cores <- ifelse(mc.cores.old>nbCores, nbCores, mc.cores.old)
        }
    } else {
        mc.cores <- 1
    }
    return(mc.cores)
}


mclapply2 <- function (X, FUN, ..., mc.preschedule = TRUE, mc.set.seed = TRUE, mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L), mc.cleanup = TRUE, mc.allow.recursive = TRUE) {
    if (Sys.info()[["sysname"]] != "Linux") {
        mc.cores <- 1
    } else {
        mc.cores <- min(detectCores(), mc.cores)
    }
    return(mclapply(X = X, FUN = FUN, ...,
                    mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed, mc.silent = mc.silent,
                    mc.cores = maxCores(mc.cores), mc.cleanup = mc.cleanup, mc.allow.recursive = mc.allow.recursive))
}


.writeSignal <- function (pattern, snpInfoDir, signalFile) {
    tmpDir <- gsub("\\\\", "/", tempdir())
    if (length(unlist(strsplit(readLines(signalFile, n = 1), split = "\t")))>1) {
        signal <- read.delim(file = signalFile, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                        colClasses = c("character", "numeric"), na.string = c("NA", ""),
                        check.names = FALSE, strip.white = TRUE, col.names =  c("SNP", "PVALUE"))
    } else {
        if (length(unlist(strsplit(readLines(signalFile, n = 1), split = " ")))>1) {
            signal <- read.delim(file = signalFile, header = TRUE, sep = " ", stringsAsFactors = FALSE,
                        colClasses = c("character", "numeric"), na.string = c("NA", ""),
                        check.names = FALSE, strip.white = TRUE, col.names =  c("SNP", "PVALUE"))
        } else {
            stop('[Enrichment:initFiles] only " " and "\t" are allowed as columns separator in Signal file.', call. = FALSE)
        }
    }

    chrom.bim <- read.delim(file = paste0(snpInfoDir, pattern, ".bim"), header = FALSE, stringsAsFactors = FALSE,
                            colClasses = c("numeric", "character", "NULL", "NULL", "NULL", "NULL"), na.string = c("NA", ""),
                            check.names = FALSE, strip.white = TRUE, col.names =  c("CHR", "SNP", "", "", "", ""))
    signal <- merge(signal, chrom.bim, by = "SNP")[, c(3, 1, 2)]
    eval(parse(text = paste0('write.table(signal, file = "', tmpDir, '/snpEnrichment/', pattern,
                                '.signal", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")')))
    return(invisible())
}


.readSignal <- function (pattern) {
    tmpDir <- gsub("\\\\", "/", tempdir())
    signal <- read.delim(file = paste0(tmpDir, "/snpEnrichment/", pattern, ".signal"), header = TRUE, stringsAsFactors = FALSE,
                            colClasses = c("NULL", "character", "numeric"), na.string = c("NA", ""),
                            check.names = FALSE, strip.white = TRUE, col.names =  c("", "SNP", "PVALUE"))
    return(signal)
}


.writeFreq <- function (pattern, snpInfoDir) {
    tmpDir <- gsub("\\\\", "/", tempdir())
    IN <- paste0(snpInfoDir, pattern)
    OUT <- paste0(tmpDir, "/snpEnrichment/", pattern)
    plinkData <- read.plink(bed = paste0(IN, ".bed"), bim = paste0(IN, ".bim"), fam = paste0(IN, ".fam"), select.snps = .readSignal(pattern)[, "SNP"])
    plinkFreq <- col.summary(plinkData$genotypes)
    plinkFreq <- cbind(snp.name = rownames(plinkFreq), MAF = plinkFreq[, "MAF"])
    plinkRes <- merge(plinkData$map, plinkFreq, by = "snp.name")
    plinkRes <- plinkRes[, c("chromosome", "snp.name", "position", "MAF")]
    plinkRes[, "MAF"] <- as.numeric(as.character(plinkRes[, "MAF"]))
    colnames(plinkRes) <- c("CHR", "SNP", "POS", "MAF")
    write.table(plinkRes, paste0(OUT, ".all"), row.names = FALSE, sep = "\t")
    return(invisible())
}


.checkSnpInfoDir <- function (snpInfoDir) {
    snpInfoDir <- .checkFilePath(snpInfoDir)
    if (length(list.files(snpInfoDir, pattern = "*.bim"))!=22) {
         stop(paste0("[Enrichment:initFiles] Only ", length(list.files(snpInfoDir, pattern = "*.bim")), " 'bim' files found when 22 is needed."), call. = FALSE)
    } else {}
    if (length(list.files(snpInfoDir, pattern = "*.bed"))!=22) {
         stop(paste0("[Enrichment:initFiles] Only ", length(list.files(snpInfoDir, pattern = "*.bed")), " 'bed' files found when 22 is needed."), call. = FALSE)
    } else {}
    if (length(list.files(snpInfoDir, pattern = "*.fam"))!=22) {
         stop(paste0("[Enrichment:initFiles] Only ", length(list.files(snpInfoDir, pattern = "*.fam")), " 'fam' files found when 22 is needed."), call. = FALSE)
    } else {}
    return(invisible())
}


.checkSignalFile <- function (signalFile) {
    if (!file.exists(signalFile)) {
         stop(paste0("[Enrichment:initFiles] ", signalFile, " doesn't exist."), call. = FALSE)
    } else {}
    return(invisible())
}


.checkSnpListDir <- function (snpListDir, pattern) {
    snpListDir <- .checkFilePath(snpListDir)
    if (length(list.files(snpListDir, pattern = pattern))==0) {
         stop(paste0("[Enrichment:readEnrichment] No snp list file found when at least one is needed."), call. = FALSE)
    } else {}
    return(invisible())
}


initFiles <- function (pattern = "Chrom", snpInfoDir, signalFile, mc.cores = 1) {
    if (missing(snpInfoDir) | missing(signalFile)) {
        stop('[Enrichment:initFiles] argument(s) missing.', call. = FALSE)
    } else {}
    snpInfoDir <- .checkFilePath(snpInfoDir)
    FILES <- list.files(snpInfoDir, pattern = ".bim")
    .checkSnpInfoDir(snpInfoDir)
    .checkSignalFile(signalFile)
    tmpDir <- gsub("\\\\", "/", tempdir())
    dir.create(paste0(tmpDir, "/snpEnrichment/"), showWarnings = FALSE)
    cat("All files are ready for chromosome:\n  ")
    resParallel <- mclapply2(X = seq(22), mc.cores = min(22, mc.cores), FUN = function (iChr) {
        # newPattern <- unlist(strsplit(grep(paste0(pattern, iChr, "[^0-9]*.bim"), FILES, value = TRUE), ".bim"))[1]
        newPattern <- gsub(".bim", "", grep(paste0(pattern, iChr, "[^0-9]"), FILES, value = TRUE))
        err1 <- try(.writeSignal(pattern = newPattern, snpInfoDir = snpInfoDir, signalFile = signalFile), silent = TRUE)
        err2 <- try(.writeFreq(pattern = newPattern, snpInfoDir = snpInfoDir), silent = TRUE)
        cat(iChr, " ", sep = "")
        if (class(err1)=="try-error" | class(err2)=="try-error") {
            return(invisible("ERROR"))
        } else {
            return(invisible())
        }
    })
    if (any(unlist(resParallel)=="ERROR")) {
        stop("[Enrichment:initFiles] initialize files failed.", call. = FALSE)
    } else {}
    cat("\n\n")
    return(invisible())
}


writeLD <- function (pattern = "Chrom", snpInfoDir, signalFile, ldDir = NULL, ldThresh = 0.8, depth = 1000, mc.cores = 1) {
    if (missing(pattern) | missing(snpInfoDir) | missing(signalFile) | missing(ldThresh)) {
        stop('[Enrichment:writeLD] argument(s) missing.', call. = FALSE)
    } else {}
    tmpDir <- gsub("\\\\", "/", tempdir())
    dir.create(paste0(tmpDir, "/snpEnrichment/"), showWarnings = FALSE)
    snpInfoDir <- .checkFilePath(snpInfoDir)
    .checkSnpInfoDir(snpInfoDir)
    .checkSignalFile(signalFile)
    if (missing(ldDir) | is.null(ldDir)) {
        ldDir <- paste0(tmpDir, "/snpEnrichment/")
    } else {
        ldDir <- .checkFilePath(ldDir)
    }
    FILES <- list.files(snpInfoDir, pattern = ".bim")
    cat("Compute LD for chromosome:\n  ")
    resParallel <- mclapply2(X = seq(22), mc.cores = min(22, mc.cores), FUN = function (iChr) {
        # newPattern <- unlist(strsplit(grep(paste0(pattern, iChr, "[^0-9]*.bim"), FILES, value = TRUE), ".bim"))[1]
        newPattern <- gsub(".bim", "", grep(paste0(pattern, iChr, "[^0-9]"), FILES, value = TRUE))
        # isThereSignals <- grep(".signal", list.files(paste0(tmpDir, "/snpEnrichment/"), full.names = TRUE), value = TRUE)
        isThereSignals <- list.files(paste0(tmpDir, "/snpEnrichment/"), full.names = TRUE, pattern = ".signal")
        if (length(isThereSignals) != 22) {
            .writeSignal(pattern = newPattern, snpInfoDir, signalFile)
        } else {}
        IN <- paste0(snpInfoDir, newPattern)
        OUT <- paste0(tmpDir, "/snpEnrichment/", newPattern)
        plinkData <- read.plink(bed = paste0(IN, ".bed"), bim = paste0(IN, ".bim"), fam = paste0(IN, ".fam"), select.snps = .readSignal(newPattern)[, "SNP"])
        ldData <- ld(x = plinkData$genotypes, depth = min(ncol(plinkData$genotypes)-1, depth), stats = "R.squared")
        if (any(isNA <- is.na(ldData))) {
            ldData <- replace(ldData, grep(TRUE, isNA), 0)
        } else {}
        ldData <- apply(ldData, 1, function (li) { which(li>ldThresh) })
        resLD <- data.frame(matrix(unlist(strsplit(names(unlist(ldData)), "\\.")), ncol = 2, byrow = TRUE, dimnames = list(NULL, c("SNP_A", "SNP_B"))), stringsAsFactors = FALSE)
        resLD <- resLD[which(resLD[, 1]!=resLD[, 2]), ]

        ### Check LD distance
        # snpNames <- grep("[sS][nN][pP]", colnames(plinkData$map), value = TRUE)
        # posNames <- grep("[pP][oO][sS].*", colnames(plinkData$map), value = TRUE)
        # tmpPlink <- plinkData$map[which(plinkData$map[, snpNames]%in%unique(unlist(resLD[, c("SNP_A", "SNP_B")]))), c(snpNames, posNames)]
        # resLD[, "Dist"] <- apply(resLD[, c("SNP_A", "SNP_B")], 1, function (li) {
            # abs(tmpPlink[tmpPlink[, snpNames]%in%li["SNP_A"], posNames] - tmpPlink[tmpPlink[, snpNames]%in%li["SNP_B"], posNames])
        # })
        # table(resLD[, "Dist"]<1000000)

        write.table(resLD[, c("SNP_A", "SNP_B")], file = paste0(ldDir, newPattern, ".ld"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        cat(iChr, " ", sep = "")
    })
    cat("\n\n")
    return(invisible())
}


.readSNP <- function (pattern, snpListDir) {
    # snpListFile <- grep(pattern, sort(list.files(snpListDir, full.names = TRUE)), value = TRUE)[1]
    snpListFile <- list.files(gsub("/$", "", snpListDir), pattern = paste0(pattern, "[^0-9]"), full.names = TRUE)
    if (is.na(snpListFile)) {
        snpList <- data.frame()
    } else {
        snpList <- read.delim(file = snpListFile, header = FALSE, stringsAsFactors = FALSE,
                                na.string = c("NA", ""), check.names = FALSE, strip.white = TRUE)
        switch(EXPR = as.character(ncol(snpList)),
            "1" = {colnames(snpList) <- "SNP"},
            "2" = {colnames(snpList) <- c("CHR", "SNP")},
            "3" = {colnames(snpList) <- c("CHR", "SNP", "TRANSCRIPT")},
            colnames(snpList)[1:3] <- c("CHR", "SNP", "TRANSCRIPT")
        )
        if (ncol(snpList)>=3) {
            snpList <- snpList[, c("SNP", "TRANSCRIPT")]
        } else {}
        if (nrow(snpList)!=0) {
            snpList[, "eSNP"] <- 1
        } else {}
    }
    return(snpList)
}


.readTranscript <- function (transcriptFile) {
    if (all(class(try(close(file(transcriptFile)), silent = TRUE))!="try-error")) {
        transcript <- read.delim(file = transcriptFile, header = TRUE, stringsAsFactors = FALSE, na.string = c("NA", ""), check.names = FALSE, strip.white = TRUE)
        transcript <- transcript[, 1:4]
        colnames(transcript) <- c("TRANSCRIPT", "CHR", "START", "END")
    } else {
        if (class(transcriptFile) %in% c("matrix", "data.frame")) {
            transcript <- transcriptFile[, 1:4]
            colnames(transcript) <- c("TRANSCRIPT", "CHR", "START", "END")
        } else {
            stop('[Enrichment:readEnrichment] "transcriptFile" required "matrix", "data.frame" or "txt file".', call. = FALSE)
        }
    }
    return(transcript)
}


.readFreq <- function (pattern, snpInfoDir) {
    tmpDir <- gsub("\\\\", "/", tempdir())
    freq <- read.delim(file = paste0(tmpDir, "/snpEnrichment/", pattern, ".all"), header = TRUE, stringsAsFactors = FALSE,
                        colClasses = c("numeric", "character", "integer", "numeric"), na.string = c("NA", ""),
                        check.names = FALSE, strip.white = TRUE, col.names =  c("CHR", "SNP", "POS", "MAF"))
    return(freq)
}


.readLD <- function (pattern, snpInfoDir, ldDir) {
    tmpDir <- gsub("\\\\", "/", tempdir())
    if (missing(ldDir) | is.null(ldDir)) {
        IN <- paste0(tmpDir, "/snpEnrichment/", pattern, ".ld")
    } else {
        IN <- paste0(.checkFilePath(ldDir), pattern, ".ld")
    }
    nbCol <- as.character(ncol(read.table(text = readLines(con = IN, n = 1), stringsAsFactors = FALSE)))
    switch(EXPR = nbCol,
        "2" = {ldData <- read.delim(file = IN, header = TRUE, stringsAsFactors = FALSE, colClasses = c("character", "character"),
                                    na.string = c("NA", ""), check.names = FALSE, strip.white = TRUE, sep = "\t")
        },
        "7" = {ldData <- read.delim(file = IN, header = TRUE, stringsAsFactors = FALSE, colClasses = c("NULL", "NULL", "character", "NULL", "NULL", "character", "NULL"),
                                    na.string = c("NA", ""), check.names = FALSE, strip.white = TRUE, sep = "")
        },
        stop(paste0('[Enrichment:readEnrichment] "', pattern, '.ld" structure must be a matrix file with columns:\n       c("SNP_A", "SNP_B") or c("CHR_A", "BP_A", "SNP_A", "CHR_B", "BP_B", "SNP_B", "R2")".'), call. = FALSE)
    )
    return(ldData)
}


.checkTranscript <- function (data, transcriptFile, distThresh) {
    if (any(transcriptFile != FALSE)) {
        transcript <- .readTranscript(transcriptFile = transcriptFile)
        transcriptCHR <- na.exclude(transcript[transcript[, "CHR"] == unique(data[, "CHR"]), c("START", "END")])

        cisFunc <- function (line, distThresh, dataTranscript) {
            position <- as.numeric(line[4])
            CIS <- any(position > (dataTranscript[, "START"] - distThresh) & (position < (dataTranscript[, "END"] + distThresh)))
            return(c(line, CIS = as.numeric(CIS)))
        }

        temp <- unique(as.data.frame(t(apply(data, 1, cisFunc, distThresh = distThresh*10^3, dataTranscript = transcriptCHR)), stringsAsFactors = FALSE))
        temp <- temp[temp[, "CIS"] == 1, -grep("CIS", colnames(temp))]
    } else {
        temp <- unique(data)
    }
    return(temp)
}


.readFiles <- function (pattern, snpInfoDir, snpListDir, distThresh) {
    fullPattern <- gsub(".bim", "", grep(paste0(pattern, "[^0-9]"), list.files(snpInfoDir, pattern = ".bim"), value = TRUE))
    signalPattern <- gsub(".signal", "", grep(paste0(pattern, "[^0-9]"), list.files(paste0(gsub("\\\\", "/", tempdir()), "/snpEnrichment/"), pattern = ".signal"), value = TRUE))
    eSNP <- .readSNP(pattern = fullPattern, snpListDir = snpListDir)
    signal <- .readSignal(pattern = signalPattern)
    plinkData <- .readFreq(pattern = fullPattern, snpInfoDir = snpInfoDir)

    signalPlink <- merge(signal, plinkData, by = "SNP")
    if (nrow(eSNP)!=0) {
        eSNPunique <- eSNP[!duplicated(eSNP[, "SNP"]), ]
        snpSignal <- merge(signalPlink, eSNPunique[, c("SNP", "eSNP")], by = "SNP", all.x = TRUE)
        snpSignal[, "eSNP"][is.na(snpSignal[, "eSNP"])] <- 0
        snpLoss <- c(
            length(eSNP[, "SNP"]),
            length(unique(eSNP[, "SNP"])),
            length(unique(intersect(eSNP[, "SNP"], signal[!is.na(signal[, "PVALUE"]), "SNP"])))
        )
    } else {
        snpSignal <- signalPlink
        snpSignal[, "eSNP"] <- 0
        snpLoss <- c(0, 0, 0)
    }

    temp <- unique(snpSignal[, c("SNP", "PVALUE", "CHR", "POS", "MAF", "eSNP")])
    data <- transform(temp, SNP = as.character(temp$SNP), PVALUE = as.numeric(temp$PVALUE), CHR = as.numeric(temp$CHR), POS = as.integer(temp$POS), MAF = as.numeric(temp$MAF), eSNP = as.numeric(temp$eSNP))
    data[, "xSNP"] <- 0
    if (any(duplicated(data[, "SNP"]))) {
        cat(data[, "SNP"][duplicated(data[, "SNP"])], "\n")
        stop('[Enrichment:readEnrichment] Duplicated SNPs in Signal.', call. = FALSE)
    } else {
        rownames(data) <- data[, "SNP"]
    }
    return(list(data = data, snpLoss = snpLoss))
}


.splitByChrom <- function (pattern, snpListFile, directory) {
    if (missing(snpListFile)) {
        stop('[Enrichment:readEnrichment] argument(s) missing.', call. = FALSE)
    } else {}
    snpList <- read.delim(file = snpListFile, header = FALSE, stringsAsFactors = FALSE,
                            na.string = c("NA", ""), check.names = FALSE, strip.white = TRUE)
    switch(EXPR = as.character(ncol(snpList)),
        "1" = {stop('[Enrichment:readEnrichment] at least two columns are needed: "Chromosome" and "rs name".', call. = FALSE)},
        "2" = {colnames(snpList) <- c("CHR", "SNP")},
        "3" = {colnames(snpList) <- c("CHR", "SNP", "TRANSCRIPT")},
        colnames(snpList)[1:3] <- c("CHR", "SNP", "TRANSCRIPT")
    )
    if (ncol(snpList)>3) {
        snpList <- snpList[, c("CHR", "SNP", "TRANSCRIPT")]
    } else {}
    filePathDetails <- unlist(strsplit(snpListFile, "/"))
    if (is.null(directory)) {
        filePath <- paste0(c(filePathDetails[-length(filePathDetails)], ""), collapse = "/")
    } else {
        dir.create(paste0(directory, "snpList/"), showWarnings = FALSE)
        filePath <- paste0(directory, "snpList/")
    }
    by(snpList, snpList[, "CHR"], function (snpListChr) {
        write.table(snpListChr, file = paste0(filePath, pattern, unique(snpListChr[, "CHR"]), "-", filePathDetails[length(filePathDetails)]), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
    })
    # cat(paste0('*** File "', filePathDetails[length(filePathDetails)], '" has been splitted into ', length(unique(snpList[, "CHR"])), ' chromosomes. ***\n\n'))
    return(filePath)
}


readEnrichment <- function (pattern = "Chrom", signalFile, transcriptFile = FALSE, snpListDir, snpInfoDir, distThresh = 1000, sigThresh = 0.05, LD = FALSE, ldDir = NULL, mc.cores = 1) {
    cat("############# Read Enrichment ##############\n")
    if (missing(signalFile) | missing(snpListDir) | missing(snpInfoDir)) {
        stop('[Enrichment:readEnrichment] argument(s) missing.', call. = FALSE)
    } else {}
    tmpDir <- paste0(gsub("\\\\", "/", tempdir()), "/snpEnrichment/")
    if (length(list.files(tmpDir, pattern = "\\.signal")) != 22 & length(list.files(tmpDir, pattern = "\\.all")) != 22) {
        stop('[Enrichment:readEnrichment] "initFiles" must be run before.', call. = FALSE)
    } else {}

    if (is.null(ldDir)) {
        if (LD & length(list.files(tmpDir, pattern = "\\.ld")) != 22) {
            stop('[Enrichment:readEnrichment] "writeLD" must be run before. (Or LD computation by PLINK.)', call. = FALSE)
        } else {}
    } else {
        if (LD & length(list.files(ldDir, pattern = "\\.ld")) != 22) {
            stop('[Enrichment:readEnrichment] "writeLD" must be run before. (Or LD computation by PLINK.)', call. = FALSE)
        } else {}
    }

    snpInfoDir <- .checkFilePath(snpInfoDir)
    .checkSnpInfoDir(snpInfoDir)
    .checkSignalFile(signalFile)

    if (length(list.files(snpListDir)) == 0) {
        snpListDir <- .splitByChrom(pattern = pattern, snpListFile = snpListDir, directory = tmpDir)
    } else {
        snpListDir <- .checkFilePath(snpListDir)
        .checkSnpListDir(snpListDir, pattern)
    }

    cat("  Read Chromosomes:\n    ")
    resParallel <- mclapply2(seq(22), mc.cores = min(22, mc.cores), FUN = function (iChr) {
        files <- .readFiles(pattern = paste0(pattern, iChr), snpInfoDir = snpInfoDir, snpListDir = snpListDir, distThresh = distThresh)
        if (LD) {
            # newPattern <- unlist(strsplit(grep(paste0(pattern, iChr, "[^0-9]*.bim"), list.files(snpInfoDir), value = TRUE), ".bim"))[1]
            newPattern <- gsub(".bim", "", grep(paste0(pattern, iChr, "[^0-9]"), list.files(snpInfoDir, pattern = ".bim"), value = TRUE))
            linkageData <- .readLD(pattern = newPattern, snpInfoDir = snpInfoDir, ldDir = ldDir)

            data <- files$data[order(files$data$POS), ]
            snpSignal <- data[, "SNP"]

            ldData <- linkageData[, 2]
            names(ldData) <- linkageData[, 1]
            eSNP <- data[data[, "eSNP"]==1, "SNP"]
            ldSNP <- linkageData[linkageData[, "SNP_A"] %in% eSNP, "SNP_B"]

            xSNP <- unique(c(eSNP, ldSNP))
            data[data[, "SNP"]%in%xSNP, "xSNP"] <- 1
            data <- data[!is.na(data[, "PVALUE"]), ]
            resChr <- chromosome(Data = data, LD = ldData)
        } else {
            resChr <- chromosome(Data = files$data)
        }
        snpLoss <- files$snpLoss

        if (all(transcriptFile == FALSE)) {
            resCheckData <- resChr@Data
            resCheckData <- transform(unique(resCheckData),
                SNP = as.character(resCheckData$SNP),
                PVALUE = as.numeric(resCheckData$PVALUE),
                CHR = as.numeric(resCheckData$CHR),
                POS = as.numeric(resCheckData$POS),
                MAF = as.numeric(resCheckData$MAF),
                eSNP = as.numeric(resCheckData$eSNP),
                xSNP = as.numeric(resCheckData$xSNP)
            )
            signalLoss <- c(NA, length(unique(resChr@Data[, "SNP"])), length(unique(resChr@Data[!is.na(resChr@Data[, "PVALUE"]), "SNP"])))
        } else {
            resCheckData <- .checkTranscript(data = resChr@Data, transcriptFile = transcriptFile, distThresh = distThresh)
            snpLoss <- c(snpLoss, sum(as.numeric(resCheckData[, "eSNP"])))
            resCheckData <- transform(unique(resCheckData),
                SNP = as.character(resCheckData$SNP),
                PVALUE = as.numeric(resCheckData$PVALUE),
                CHR = as.numeric(resCheckData$CHR),
                POS = as.numeric(resCheckData$POS),
                MAF = as.numeric(resCheckData$MAF),
                eSNP = as.numeric(resCheckData$eSNP),
                xSNP = as.numeric(resCheckData$xSNP)
            )
            signalLoss <- c(NA, length(unique(resChr@Data[, "SNP"])), length(unique(resChr@Data[!is.na(resChr@Data[, "PVALUE"]), "SNP"])), length(unique(resCheckData[!is.na(resCheckData[, "PVALUE"]), "SNP"])))
        }

        resChr <- chromosome(Data = resCheckData, LD = resChr@LD)

        cat(iChr, " ", sep = "")
        return(list(resChr, snpLoss, signalLoss))
    })
    names(resParallel) <- paste0("Chrom", seq(22))

    result = enrichment()
    result@Chromosomes <- lapply(resParallel, "[[", 1)
    snpLoss <- t(sapply(resParallel, "[[", 2))
    signalLoss <- rowSums(sapply(resParallel, "[[", 3))
    loss <- rbind(Signal = signalLoss, Genome = apply(snpLoss, 2, sum), snpLoss)
    colnames(loss) <- c("Rows", "Unique", "Intersect.Signal", "CIS")[seq(ncol(loss))]

    SNPs <- result["List", seq(22)]
    result@eSNP@List <- SNPs[["eSNP"]]
    result@xSNP@List <- SNPs[["xSNP"]]

    if (length(unlist(strsplit(readLines(signalFile, n = 1), split = "\t")))>1) {
        signal <- read.delim(file = signalFile, header = TRUE, sep = "\t", colClasses = c("character", "numeric"), na.string = c("NA", ""), check.names = FALSE, strip.white = TRUE, col.names =  c("SNP", "PVALUE"), stringsAsFactors = FALSE)
    } else {
        if (length(unlist(strsplit(readLines(signalFile, n = 1), split = " ")))>1) {
            signal <- read.delim(file = signalFile, header = TRUE, sep = " ", colClasses = c("character", "numeric"), na.string = c("NA", ""), check.names = FALSE, strip.white = TRUE, col.names =  c("SNP", "PVALUE"), stringsAsFactors = FALSE)
        } else {
            stop('[Enrichment:readEnrichment] only " " and "\t" are allowed as columns separator in Signal file.', call. = FALSE)
        }
    }
    loss[1, c(1, 2)] <- c(length(signal[, "SNP"]), length(unique(signal[, "SNP"])))
    result@Loss <- as.data.frame(loss)

    sysCall <- sys.call(sys.parent())
    argsSNP <- as.list(sysCall[-1])
    formal <- as.list(names(formals(as.character(sysCall))))
    names(formal) <- formal
    if (is.null(names(argsSNP))) {
        names(argsSNP) <- names(formal)[seq_along(argsSNP)]
    } else {
        emptyNames <- which(names(argsSNP)=="")
        names(argsSNP)[emptyNames] <- names(formal)[emptyNames]
    }
    names(argsSNP)[grep("^$", names(argsSNP))] <- names(formal)[grep("^$", names(argsSNP))]
    argsSNP <- c(argsSNP, lapply(formal[!names(formal) %in% names(argsSNP)], as.name))[names(formal)]
    paramClass <- sapply(argsSNP, class)

    for (iArg in names(formal)) {
        if (iArg != "...") {
            paramPos <- grep(iArg, names(formal), fixed = TRUE)
            argTmp <- argsSNP[[paramPos]]
            classTmp <- paramClass[[paramPos]]
            switch(EXPR = classTmp,
                "character" = {formal[[iArg]] <- argTmp},
                "logical" = {formal[[iArg]] <- argTmp},
                "numeric" = {formal[[iArg]] <- argTmp},
                "integer" = {formal[[iArg]] <- argTmp},
                "NULL" = {formal[[iArg]] <- "NULL"},
                "call" = {formal[[iArg]] <- eval(argTmp)},
                "name" = {
                    resEval <- eval(argTmp)
                    switch(EXPR = class(resEval),
                        "character" = {formal[[iArg]] <- resEval},
                        "logical" = {formal[[iArg]] <- resEval},
                        "numeric" = {formal[[iArg]] <- resEval},
                        "integer" = {formal[[iArg]] <- resEval},
                        "matrix" = {formal[[iArg]] <- argTmp},
                        "data.frame" = {formal[[iArg]] <- argTmp},
                        "Enrichment" = {formal[[iArg]] <- argTmp}
                    )
                }
            )
        } else {}
    }

    # unlink(tmpDir, recursive = TRUE) # add parameter to ask if temporary directory should be cleaned
    if (is.null(ldDir)) {
        formal[["ldDir"]] <- tmpDir
    } else {
        formal[["ldDir"]] <- .checkFilePath(formal[["ldDir"]])
    }
    result@Call$readEnrichment <- formal

    result <- computeER(object = result, sigThresh = sigThresh, mc.cores = mc.cores)
    cat("\n########### Read Enrichment Done ###########\n\n")
    return(result)
}


.enrichmentRatio <- function (table) {
    tmp <- as.numeric(table)
    return((tmp[1]*tmp[4])/(tmp[2]*tmp[3]))
}


.reSample <- function (object, nSample, empiricPvalue, sigThresh, MAFpool, mc.cores) {
    nResampling <- nrow(object@eSNP@Resampling)
    data <- object["Data"]
    chrLD <- object["LD"]
    isLD <- length(chrLD) != 0
    eEnrichment <- object@eSNP@Table
    eEnrichRatio <- object@eSNP@EnrichmentRatio
    eSNPlist <- object@eSNP@List
    data[, "MAFpool"] <- NA
    data[, "MAFpool"] <- as.factor(cut(data[, "MAF"], breaks = MAFpool, labels = FALSE, include.lowest = TRUE))
    nPool <- nlevels(data[, "MAFpool"])
    eSNPlistPool <- table(data[data[, "eSNP"] == 1, "MAFpool"])

    if (nResampling == 0) {
        eSNPsample <- NULL
        nSampleMin <- min(nSample, max(1000, nSample*10/100))
    } else {
        eSNPsample <- object@eSNP@Resampling
        nSampleMin <- nSample
    }

    if (isLD) {
        if (nResampling == 0) {
            xSNPsample <- NULL
        } else {
            xSNPsample <- object@xSNP@Resampling
        }
        xEnrichment <- object@xSNP@Table
        xEnrichRatio <- object@xSNP@EnrichmentRatio
        xSNPlist <- object@xSNP@List
        xSNPlistPool <- table(data[data[, "xSNP"] == 1, "MAFpool"], deparse.level = 0)
        popSNP4Sample <- split(data[data[, "xSNP"] != 1, "SNP"], data[data[, "xSNP"] != 1, "MAFpool"])
    } else {
        popSNP4Sample <- split(data[data[, "eSNP"] != 1, "SNP"], data[data[, "eSNP"] != 1, "MAFpool"])
    }
    assoc <- factor(data[, "PVALUE"] < sigThresh, levels = c(FALSE, TRUE))
    SNPlist <- data[, "SNP"]
    rm(data)

    cat(0, ".. ", sep = "")
    resParallel <- mclapply2(X = seq(nSampleMin), mc.cores = mc.cores, FUN = function (i) {
        eSNPlistRandom <- unlist(sapply(seq(nPool), function (g) {sample(popSNP4Sample[[g]], min(eSNPlistPool[g], length(popSNP4Sample[[g]])))}))
        eSNPenrichStats <- table(assoc, factor(SNPlist%in%eSNPlistRandom, levels = c(FALSE, TRUE)), deparse.level = 0)
        if (isLD) {
            xSNPlistRandom <- intersect(SNPlist, unique(chrLD[which(names(chrLD)%in%eSNPlistRandom)]))
            xSNPenrichStats <- table(assoc, factor(SNPlist%in%xSNPlistRandom, levels = c(FALSE, TRUE)), deparse.level = 0)
            xTmp <- c(c(xSNPenrichStats), .enrichmentRatio(xSNPenrichStats))
        } else {
            xTmp <- NULL
        }
        eTmp <- c(c(eSNPenrichStats), .enrichmentRatio(eSNPenrichStats))
        return(list(eSNP = eTmp, xSNP = xTmp))
    })
    eSNPsample <- rbind(eSNPsample, do.call("rbind", lapply(resParallel, function (l) {l$eSNP})))
    eMeanEnrichRatio <- sum(eSNPsample[, 5])/length(eSNPsample[, 5])
    Ze <- (eEnrichRatio-eMeanEnrichRatio)/sqrt((sum((eSNPsample[, 5]-eMeanEnrichRatio)^2))/(length(eSNPsample[, 5])-1))
    if (isLD) {
        xSNPsample <- rbind(xSNPsample, do.call("rbind", lapply(resParallel, function (l) {l$xSNP})))
        xMeanEnrichRatio <- sum(xSNPsample[, 5])/length(xSNPsample[, 5])
        Zx <- (xEnrichRatio-xMeanEnrichRatio)/sqrt((sum((xSNPsample[, 5]-xMeanEnrichRatio)^2))/(length(xSNPsample[, 5])-1))
    } else {
        Zx <- 1
    }
    iSample <- nSampleMin
    cat(iSample, ".. ", sep = "")
    catStep <- (nSample/10)
    rm(resParallel)

    while (iSample<nSample) {
        resParallel <- mclapply2(X = seq(nSampleMin), mc.cores = mc.cores, FUN = function (i) {
            eSNPlistRandom <- unlist(sapply(seq(nPool), function (g) {sample(popSNP4Sample[[g]], min(eSNPlistPool[g], length(popSNP4Sample[[g]])))}))
            eSNPenrichStats <- table(assoc, factor(SNPlist%in%eSNPlistRandom, levels = c(FALSE, TRUE)), deparse.level = 0)
            if (isLD) {
                xSNPlistRandom <- intersect(SNPlist, unique(chrLD[which(names(chrLD)%in%eSNPlistRandom)]))
                xSNPenrichStats <- table(assoc, factor(SNPlist%in%xSNPlistRandom, levels = c(FALSE, TRUE)), deparse.level = 0)
                xTmp <- c(c(xSNPenrichStats), .enrichmentRatio(xSNPenrichStats))
            } else {
                xTmp <- NULL
            }
            eTmp <- c(c(eSNPenrichStats), .enrichmentRatio(eSNPenrichStats))
            return(list(eSNP = eTmp, xSNP = xTmp))
        })
        eSNPsample <- rbind(eSNPsample, do.call("rbind", lapply(resParallel, function (l) {l$eSNP})))
        eMeanEnrichRatio <- sum(eSNPsample[, 5])/length(eSNPsample[, 5])
        Ze <- (eEnrichRatio-eMeanEnrichRatio)/sqrt((sum((eSNPsample[, 5]-eMeanEnrichRatio)^2))/(length(eSNPsample[, 5])-1))
        if (isLD) {
            xSNPsample <- rbind(xSNPsample, do.call("rbind", lapply(resParallel, function (l) {l$xSNP})))
            xMeanEnrichRatio <- sum(xSNPsample[, 5])/length(xSNPsample[, 5])
            Zx <- (xEnrichRatio-xMeanEnrichRatio)/sqrt((sum((xSNPsample[, 5]-xMeanEnrichRatio)^2))/(length(xSNPsample[, 5])-1))
        } else {
            Zx <- 1
        }
        iSample <- iSample + nSampleMin
        if ((iSample-(iSample%/%catStep*catStep)) == 0) {cat(iSample, ".. ", sep = "")}
        rm(resParallel)
    }

    if (empiricPvalue) {
        whichPvalue <- 2
    } else {
        whichPvalue <- 1
    }

    empiricPvalue.eSNP <- sum(eEnrichRatio<eSNPsample[, 5])/length(eSNPsample[, 5])
    statisticPvalue.eSNP <- pnorm(Ze, lower.tail = FALSE)
    object@eSNP <- enrichSNP(List = eSNPlist, Table = eEnrichment, EnrichmentRatio = eEnrichRatio, Z = Ze, PValue = c(Distribution = statisticPvalue.eSNP, Empirical = empiricPvalue.eSNP)[whichPvalue], Resampling = eSNPsample)
    if (isLD) {
        empiricPvalue.xSNP <- sum(xEnrichRatio<xSNPsample[, 5])/length(xSNPsample[, 5])
        statisticPvalue.xSNP <- pnorm(Zx, lower.tail = FALSE)
        object@xSNP <- enrichSNP(List = xSNPlist, Table = xEnrichment, EnrichmentRatio = xEnrichRatio, Z = Zx, PValue = c(Distribution = statisticPvalue.xSNP, Empirical = empiricPvalue.xSNP)[whichPvalue], Resampling = xSNPsample)
    } else {
        object@xSNP <- enrichSNP()
    }
    return(object)
}


.compareEnrich <- function (object1, object2, nSample, empiricPvalue, sigThresh, MAFpool, mc.cores) {
    DATA <- object1["Data"]
    chrLD <- object1["LD"]
    isLD <- length(chrLD) != 0
    eEnrichment <- object1@eSNP@Table
    eEnrichRatio <- object1@eSNP@EnrichmentRatio
    eSNPlist <- object1@eSNP@List
    eList <- union(object1@eSNP@List, object2@eSNP@List)
    if (isLD) {
        xList <- union(object1@xSNP@List, object2@xSNP@List)
        data <- DATA[DATA[, "SNP"] %in% union(eList, xList), ]
    } else {
        data <- DATA[DATA[, "SNP"] %in% eList, ]
    }

    data[, "MAFpool"] <- NA
    data[, "MAFpool"] <- as.factor(cut(data[, "MAF"], breaks = MAFpool, labels = FALSE, include.lowest = TRUE))
    nPool <- nlevels(data$MAFpool)
    eSNPlistPool <- table(data[eSNPlist, "MAFpool"])

    eSNPsample <- NULL
    nSampleMin <- min(nSample, max(1000, nSample*10/100))
    if (isLD) {
        xEnrichment <- object1@xSNP@Table
        xEnrichRatio <- object1@xSNP@EnrichmentRatio
        xSNPlist <- object1@xSNP@List
        xSNPlistPool <- table(data[xSNPlist, "MAFpool"])
        xSNPsample <- NULL
    } else {}

    popSNP4Sample <- split(eList, data[eList, "MAFpool"])
    assoc_eSNP <- factor(data[eList, "PVALUE"] < sigThresh, levels = c(FALSE, TRUE))
    if (isLD) {
        assoc_xSNP <- factor(data[xList, "PVALUE"] < sigThresh, levels = c(FALSE, TRUE))
    } else { }
    rm(data)

    cat("0.. ")
    resParallel <- mclapply2(X = seq(nSampleMin), mc.cores = mc.cores, FUN = function (i) {
        eSNPlistRandom <- unlist(sapply(seq(nPool), function (g) {sample(popSNP4Sample[[g]], min(eSNPlistPool[g], length(popSNP4Sample[[g]])))}))
        eSNPenrichStats <- table(assoc_eSNP, factor(eList%in%eSNPlistRandom, levels = c(FALSE, TRUE)))
        eTmp <- c(c(eSNPenrichStats), .enrichmentRatio(eSNPenrichStats))
        if (isLD) {
            xSNPlistRandom <- intersect(xList, unique(chrLD[which(names(chrLD)%in%eSNPlistRandom)]))
            xSNPenrichStats <- table(assoc_xSNP, factor(xList%in%xSNPlistRandom, levels = c(FALSE, TRUE)))
            xTmp <- c(c(xSNPenrichStats), .enrichmentRatio(xSNPenrichStats))
        } else {
            xTmp <- NULL
        }
        return(list(eSNP = eTmp, xSNP = xTmp))
    })
    eSNPsample <- rbind(eSNPsample, do.call("rbind", lapply(resParallel, function (l) {l$eSNP})))

    eSNPanyDup <- eSNPsample[!duplicated(eSNPsample), ]
    if (is.matrix(eSNPanyDup)) {
        if (length(which(is.na(eSNPsample[, 5]))) != nrow(eSNPsample) & length(which(is.na(eSNPsample[, 5]))) > 0) {
            eSNPsample[is.na(eSNPsample[, 5]), 5] <- rep(mean(eSNPsample[(!is.na(eSNPsample[, 5]) & is.finite(eSNPsample[, 5])), 5]), length(eSNPsample[is.na(eSNPsample[, 5]), 5]))
        } else {}
        if (length(which(is.infinite(eSNPsample[, 5]))) != nrow(eSNPsample) & length(which(is.infinite(eSNPsample[, 5]))) > 0) {
            eSNPsample[is.infinite(eSNPsample[, 5]), 5] <- rep(max(eSNPsample[!is.infinite(eSNPsample[, 5]), 5])*1.05, length(eSNPsample[is.infinite(eSNPsample[, 5]), 5]))
        } else {}
        eMeanEnrichRatio <- sum(eSNPsample[, 5])/length(eSNPsample[, 5])
        Ze <- (eEnrichRatio-eMeanEnrichRatio)/sqrt((sum((eSNPsample[, 5]-eMeanEnrichRatio)^2))/(length(eSNPsample[, 5])-1))
    } else {
        Ze <- as.numeric(NA)
    }

    if (isLD) {
        xSNPsample <- rbind(xSNPsample, do.call("rbind", lapply(resParallel, function (l) {l$xSNP})))
        xSNPanyDup <- xSNPsample[!duplicated(xSNPsample), ]
        if (is.matrix(xSNPanyDup)) {
            if (length(which(is.na(xSNPsample[, 5]))) != nrow(xSNPsample) & length(which(is.na(xSNPsample[, 5]))) > 0) {
                xSNPsample[is.na(xSNPsample[, 5]), 5] <- rep(mean(xSNPsample[(!is.na(xSNPsample[, 5]) & is.finite(xSNPsample[, 5])), 5]), length(xSNPsample[is.na(xSNPsample[, 5]), 5]))
            } else {}
            if (length(which(is.infinite(xSNPsample[, 5]))) != nrow(xSNPsample) & length(which(is.infinite(xSNPsample[, 5]))) > 0) {
                xSNPsample[is.infinite(xSNPsample[, 5]), 5] <- rep(max(xSNPsample[!is.infinite(xSNPsample[, 5]), 5])*1.05, length(xSNPsample[is.infinite(xSNPsample[, 5]), 5]))
            } else {}
            xMeanEnrichRatio <- sum(xSNPsample[, 5])/length(xSNPsample[, 5])
            Zx <- (xEnrichRatio-xMeanEnrichRatio)/sqrt((sum((xSNPsample[, 5]-xMeanEnrichRatio)^2))/(length(xSNPsample[, 5])-1))
        } else {
            Zx <- as.numeric(NA)
        }
    } else {
        Zx <- as.numeric(NA)
    }
    iSample <- nSampleMin
    cat(iSample, ".. ", sep = "")
    catStep <- (nSample/10)
    rm(resParallel)

    while (iSample<nSample) {
        resParallel <- mclapply2(X = seq(nSampleMin), mc.cores = mc.cores, FUN = function (i) {
            eSNPlistRandom <- unlist(sapply(seq(nPool), function (g) {sample(popSNP4Sample[[g]], min(eSNPlistPool[g], length(popSNP4Sample[[g]])))}))
            eSNPenrichStats <- table(assoc_eSNP, factor(eList%in%eSNPlistRandom, levels = c(FALSE, TRUE)))
            eTmp <- c(c(eSNPenrichStats), .enrichmentRatio(eSNPenrichStats))
            if (isLD) {
                xSNPlistRandom <- intersect(xList, unique(chrLD[which(names(chrLD)%in%eSNPlistRandom)]))
                xSNPenrichStats <- table(assoc_xSNP, factor(xList%in%xSNPlistRandom, levels = c(FALSE, TRUE)))
                xTmp <- c(c(xSNPenrichStats), .enrichmentRatio(xSNPenrichStats))
            } else {
                xTmp <- NULL
            }
            return(list(eSNP = eTmp, xSNP = xTmp))
        })
        eSNPsample <- rbind(eSNPsample, do.call("rbind", lapply(resParallel, function (l) {l$eSNP})))

        eSNPanyDup <- eSNPsample[!duplicated(eSNPsample), ]
        if (is.matrix(eSNPanyDup)) {
            if (length(which(is.na(eSNPsample[, 5]))) != nrow(eSNPsample) & length(which(is.na(eSNPsample[, 5]))) > 0) {
                eSNPsample[is.na(eSNPsample[, 5]), 5] <- rep(mean(eSNPsample[(!is.na(eSNPsample[, 5]) & is.finite(eSNPsample[, 5])), 5]), length(eSNPsample[is.na(eSNPsample[, 5]), 5]))
            } else {}
            if (length(which(is.infinite(eSNPsample[, 5]))) != nrow(eSNPsample) & length(which(is.infinite(eSNPsample[, 5]))) > 0) {
                eSNPsample[is.infinite(eSNPsample[, 5]), 5] <- rep(max(eSNPsample[is.finite(eSNPsample[, 5]), 5])*1.05, length(eSNPsample[is.infinite(eSNPsample[, 5]), 5]))
            } else {}
            eMeanEnrichRatio <- sum(eSNPsample[, 5])/length(eSNPsample[, 5])
            Ze <- (eEnrichRatio-eMeanEnrichRatio)/sqrt((sum((eSNPsample[, 5]-eMeanEnrichRatio)^2))/(length(eSNPsample[, 5])-1))
        } else {
            Ze <- as.numeric(NA)
        }

        if (isLD) {
            xSNPsample <- rbind(xSNPsample, do.call("rbind", lapply(resParallel, function (l) {l$xSNP})))
            xSNPanyDup <- xSNPsample[!duplicated(xSNPsample), ]
            if (is.matrix(xSNPanyDup)) {
                if (length(which(is.na(xSNPsample[, 5]))) != nrow(xSNPsample) & length(which(is.na(xSNPsample[, 5]))) > 0) {
                    xSNPsample[is.na(xSNPsample[, 5]), 5] <- rep(mean(xSNPsample[(!is.na(xSNPsample[, 5]) & is.finite(xSNPsample[, 5])), 5]), length(xSNPsample[is.na(xSNPsample[, 5]), 5]))
                } else {}
                if (length(which(is.infinite(xSNPsample[, 5]))) != nrow(xSNPsample) & length(which(is.infinite(xSNPsample[, 5]))) > 0) {
                    xSNPsample[is.infinite(xSNPsample[, 5]), 5] <- rep(max(xSNPsample[is.finite(xSNPsample[, 5]), 5])*1.05, length(xSNPsample[is.infinite(xSNPsample[, 5]), 5]))
                } else {}
                xMeanEnrichRatio <- sum(xSNPsample[, 5])/length(xSNPsample[, 5])
                Zx <- (xEnrichRatio-xMeanEnrichRatio)/sqrt((sum((xSNPsample[, 5]-xMeanEnrichRatio)^2))/(length(xSNPsample[, 5])-1))
            } else {
                Zx <- as.numeric(NA)
            }
        } else {
            Zx <- as.numeric(NA)
        }
        iSample <- iSample + nSampleMin
        if ((iSample-(iSample%/%catStep*catStep)) == 0) {cat(iSample, ".. ", sep = "")}
        rm(resParallel)
    }

    if (empiricPvalue) {
        whichPvalue <- 2
    } else {
        whichPvalue <- 1
    }

    empiricPvalue.eSNP <- sum(eEnrichRatio<eSNPsample[, 5])/length(eSNPsample[, 5])
    statisticPvalue.eSNP <- pnorm(Ze, lower.tail = FALSE)
    object1@eSNP <- enrichSNP(List = eSNPlist, Table = eEnrichment, EnrichmentRatio = eEnrichRatio, Z = Ze, PValue = c(Distribution = statisticPvalue.eSNP, Empirical = empiricPvalue.eSNP)[whichPvalue], Resampling = eSNPsample)
    if (isLD) {
        empiricPvalue.xSNP <- sum(xEnrichRatio<xSNPsample[, 5])/length(xSNPsample[, 5])
        statisticPvalue.xSNP <- pnorm(Zx, lower.tail = FALSE)
        object1@xSNP <- enrichSNP(List = xSNPlist, Table = xEnrichment, EnrichmentRatio = xEnrichRatio, Z = Zx, PValue = c(Distribution = statisticPvalue.xSNP, Empirical = empiricPvalue.xSNP)[whichPvalue], Resampling = xSNPsample)
    } else {
        object1@xSNP <- enrichSNP()
    }
    if (is.enrichment(object1)) {
        object1@Chromosomes <- lapply(object1@Chromosomes, reset, "Z")
        object1@Chromosomes <- lapply(object1@Chromosomes, reset, "PValue")
        object1@Chromosomes <- lapply(object1@Chromosomes, reset, "Resampling")
    } else {}
    return(object1)
}
