loadXHMMdata <- function(basename, EXCLUDE_LARGE_MATRICES=c()) {
	writeLines(paste("Loading data from '", basename, "' XHMM run...", sep=""))

	dataList = list()

	GC = loadNamedVectorNoHeaderMayNotExist(paste(basename, ".locus_GC.txt", sep=""))
	dataList[["GC"]] = GC

	Complexity = loadNamedVectorNoHeaderMayNotExist(paste(basename, ".locus_complexity.txt", sep=""))
	dataList[["Repeat-masked"]] = Complexity

	RD = readNamedMatrix(paste(basename, ".RD.txt", sep=""))
	if (!("RD" %in% EXCLUDE_LARGE_MATRICES)) {
		dataList[["RD"]] = RD
	}

	TARGETS = colnames(RD)
	remove(RD)
	TARGETS_CHR_BP1_BP2 = targetsToChrBp1Bp2(TARGETS)
	TARGET_SIZES = targetsToSizes(TARGETS)

	dataList[["TARGETS"]] = TARGETS
	dataList[["TARGETS_CHR_BP1_BP2"]] = TARGETS_CHR_BP1_BP2
	dataList[["Target size"]] = TARGET_SIZES

	dataList[["CHROMOSOMES_START_WITH_CHR"]] = length(grep("^chr", TARGETS_CHR_BP1_BP2[["chr"]], perl=TRUE)) > 0

	if (!("filtered_centered_RD" %in% EXCLUDE_LARGE_MATRICES)) {
		filtered_centered_RD = readNamedMatrix(paste(basename, ".filtered_centered.RD.txt", sep=""))
		dataList[["filtered_centered_RD"]] = filtered_centered_RD
	}

	filtered_centered_RD.filtered_samples = scanVectorMayNotExist(paste(basename, ".filtered_centered.RD.txt.filtered_samples.txt", sep=""))
	dataList[["filtered_centered_RD.filtered_samples"]] = filtered_centered_RD.filtered_samples

	filtered_centered_RD.filtered_targets = scanVectorMayNotExist(paste(basename, ".filtered_centered.RD.txt.filtered_targets.txt", sep=""))
	dataList[["filtered_centered_RD.filtered_targets"]] = filtered_centered_RD.filtered_targets

	if (!("PC" %in% EXCLUDE_LARGE_MATRICES)) {
		PC = readNamedMatrix(paste(basename, ".RD_PCA.PC.txt", sep=""))
		dataList[["PC"]] = PC
	}

	PC_SD = readNamedMatrix(paste(basename, ".RD_PCA.PC_SD.txt", sep=""))
	dataList[["PC_SD"]] = PC_SD

	PC_LOADINGS = readNamedMatrix(paste(basename, ".RD_PCA.PC_LOADINGS.txt", sep=""))
	dataList[["PC_LOADINGS"]] = PC_LOADINGS

	PCA_NORMALIZE_NUM_REMOVED = paste(basename, ".PCA_normalized.txt.num_removed_PC.txt", sep="")
	NUM_PC_REMOVED_VEC = scanVectorMayNotExist(PCA_NORMALIZE_NUM_REMOVED)
	if (length(NUM_PC_REMOVED_VEC) == 1) {
		NUM_PC_REMOVED = as.numeric(NUM_PC_REMOVED_VEC[1])
	}
	else {
		PCA_NORMALIZE_OUT_FILE = paste(basename, ".PCA_normalized.txt.out", sep="")
		if (!file.exists(PCA_NORMALIZE_OUT_FILE)) {
			stop(paste("Cannot find output log file '", PCA_NORMALIZE_OUT_FILE, "' either", sep=""))
		}
		NUM_PC_REMOVED = as.numeric(system(paste("grep 'Removing first' ", PCA_NORMALIZE_OUT_FILE, " | awk '{print $3}'", sep=""), intern=TRUE))
		if (length(NUM_PC_REMOVED) != 1) {
			stop(paste("Unable to determine number of PC removed from '", PCA_NORMALIZE_OUT_FILE, "'", sep=""))
		}
	}
	dataList[["NUM_PC_REMOVED"]] = NUM_PC_REMOVED

	if (!("PCA_NORMALIZED" %in% EXCLUDE_LARGE_MATRICES)) {
		PCA_NORMALIZED = readNamedMatrix(paste(basename, ".PCA_normalized.txt", sep=""))
		dataList[["PCA_NORMALIZED"]] = PCA_NORMALIZED
	}

	if (!("PCA_NORM_Z_SCORES" %in% EXCLUDE_LARGE_MATRICES)) {
		PCA_NORM_Z_SCORES = readNamedMatrix(paste(basename, ".PCA_normalized.filtered.sample_zscores.RD.txt", sep=""))
		dataList[["PCA_NORM_Z_SCORES"]] = PCA_NORM_Z_SCORES
	}

	PCA_NORM_Z_SCORES.filtered_samples = scanVectorMayNotExist(paste(basename, ".PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt", sep=""))
	dataList[["PCA_NORM_Z_SCORES.filtered_samples"]] = PCA_NORM_Z_SCORES.filtered_samples

	PCA_NORM_Z_SCORES.filtered_targets = scanVectorMayNotExist(paste(basename, ".PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt", sep=""))
	dataList[["PCA_NORM_Z_SCORES.filtered_targets"]] = PCA_NORM_Z_SCORES.filtered_targets

	if (!("RD_SAME_FILTERED" %in% EXCLUDE_LARGE_MATRICES)) {
		RD_SAME_FILTERED = readNamedMatrix(paste(basename, ".same_filtered.RD.txt", sep=""))
		dataList[["RD_SAME_FILTERED"]] = RD_SAME_FILTERED
	}

	if (!("DIP_POST" %in% EXCLUDE_LARGE_MATRICES)) {
		DIP_POST = readNamedMatrix(paste(basename, ".posteriors.DIP.txt", sep=""))
		dataList[["DIP_POST"]] = DIP_POST
	}

	if (!("DEL_POST" %in% EXCLUDE_LARGE_MATRICES)) {
		DEL_POST = readNamedMatrix(paste(basename, ".posteriors.DEL.txt", sep=""))
		dataList[["DEL_POST"]] = DEL_POST
	}

	if (!("DUP_POST" %in% EXCLUDE_LARGE_MATRICES)) {
		DUP_POST = readNamedMatrix(paste(basename, ".posteriors.DUP.txt", sep=""))
		dataList[["DUP_POST"]] = DUP_POST
	}

	XCNV_CALLS = loadXCNVcalls(paste(basename, ".xcnv", sep=""))
	dataList[["XCNV_CALLS"]] = XCNV_CALLS

	return(dataList)
}

targetsToChrBp1Bp2 <- function(targets, startStopDelim="-") {
	chr_basePairInterval = strsplit(targets, ":", fixed=TRUE)
        chr = sapply(chr_basePairInterval, "[[", 1)
        
        bpInterval <- sapply(chr_basePairInterval, "[[", 2)
        bp1_bp2 = strsplit(bpInterval, startStopDelim, fixed=TRUE)
        bp1 = as.numeric(sapply(bp1_bp2, "[[", 1))
        bp2 = as.numeric(sapply(bp1_bp2, function (x) {if (length(x) == 2) return(x[2]) else return(x[1])}))

	return(list(chr=chr, bp1=bp1, bp2=bp2))
}

targetsToSizes <- function(targets) {
	chrBp1Bp2 = targetsToChrBp1Bp2(targets)
	sizes = chrBp1Bp2[["bp2"]] - chrBp1Bp2[["bp1"]] + 1
	names(sizes) = targets
	return(sizes)
}

# Find segments in which ANY of the vectors in valuesVectorList change value:
calcSegmentIntervals <- function(valuesVectorList) {
	firstVec = valuesVectorList[[1]]

	twoToEnd = 2:length(firstVec)
	oneToEndMinusOne = 1:(length(firstVec)-1)

	stopIntervalInds = rep(FALSE, length(firstVec)-1)
	for (i in 1:length(valuesVectorList)) {
		vec = valuesVectorList[[i]]
		if (length(vec) != length(firstVec)) {
			stop("length(vec) != length(firstVec)")
		}
		stopIntervalInds = stopIntervalInds | as.vector(vec[twoToEnd] != vec[oneToEndMinusOne])
	}
	stopIntervalInds = which(stopIntervalInds)

	startIntervalInds = c(1, stopIntervalInds + 1)
	stopIntervalInds = c(stopIntervalInds, length(firstVec))

	cbind(startIntervalInds, stopIntervalInds)
}

plotAllChromosomeValues <- function(chr, values, ylab, title, outputName, specialColors=NULL) {
	startEndChrIntervals = calcSegmentIntervals(list(chr))
        startChrinds = startEndChrIntervals[, 1]
        endChrinds = startEndChrIntervals[, 2]

	chrLabels = chr[startChrinds]
	chrColors = as.vector(matrix(nrow=length(chr)))
	for (j in 1:nrow(startEndChrIntervals)) {
		segment = startEndChrIntervals[j,1]:startEndChrIntervals[j,2]
		chrColors[segment] = j
	}
	chrColors = chrColors %% 2
	chrColors[chrColors == 1] = "black"
	chrColors[chrColors == 0] = "gray"

	colors = chrColors
	if (!is.null(specialColors)) {
		useSpecialInds = which(!is.na(specialColors))
		colors[useSpecialInds] = specialColors[useSpecialInds]
	}

	png(paste(outputName, ".png", sep=""), width=14, height=7, units="in", res=100)

	plot(values, main=title, cex.main=0.85, xaxt="n", col=colors, xlab="Exome target", ylab=ylab, pch=".", cex=3)
	axis(1, at=startChrinds, labels=chrLabels, cex.axis=0.6)

	dev.off()
}


listOfNulls <- function(n) {
        l = list()
        l[[n+1]] = "a"
        l[[n+1]] = NULL
        return(l)
}


############################
# AS IN XHMM C++ CODE:
# Instead of splitting by whitespace, use ONLY tab as delimiter (to allow for sample names with space in them):
############################

loadXCNVcalls <- function(xcnvFile) {
	return(read.table(xcnvFile, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, colClasses="character", sep="\t"))
}
