DEFAULT_MARK_SAMPLES_COLORS=c("darkorchid1", "darkorange", "cyan3")

plot_XHMM_region <- function(BASE_OUTPUT_NAME, xhmm_data, allTargsToGenes, binarySampleFeatures, SQ_THRESH, chrom, start=NULL, stop=NULL, PLOT_ONLY_PNG=TRUE, PLOT_LINE_SEGMENTS=TRUE, DEFAULT_WIDTH=14, DEFAULT_HEIGHT=7, UPDATE_FIG_HEIGHT_FOR_MARGINS=TRUE, GENE_NAME_cex=0.8, PLOT_PIPELINE_TRANSITIONS=TRUE, DIFFERENTIATE_SAMPLE_GROUPS=NULL, COLOR_DIFFERENTIAL_TARGETS=TRUE, MARK_SAMPLES=NULL, MARK_SAMPLES_COLORS=DEFAULT_MARK_SAMPLES_COLORS, MARK_SAMPLES_ONLY_IN_CALLS=FALSE, EXCLUDE_SAMPLES=NULL, MARK_INTERVALS=NULL, APPEND_REGION_NAME=TRUE) {
	chrom_start_stop = parseRegion(chrom, start, stop)
	chrom = as.character(chrom_start_stop["chrom"])
	start = as.numeric(chrom_start_stop["start"])
	stop  = as.numeric(chrom_start_stop["stop"])

	targFinalInds = regionToFinalInds(xhmm_data, chrom, start, stop)
	targs = colnames(xhmm_data[["PCA_NORM_Z_SCORES"]])[targFinalInds]

	OUTPUT_NAME = BASE_OUTPUT_NAME
	if (APPEND_REGION_NAME) {
		OUTPUT_NAME = paste(OUTPUT_NAME, ".region_", chrom, ":", format(start, scientific=FALSE), "-", format(stop, scientific=FALSE), sep="")
	}

	plot_XHMM_targets_helper(OUTPUT_NAME, xhmm_data, allTargsToGenes, binarySampleFeatures, SQ_THRESH, targs, PLOT_ONLY_PNG, PLOT_LINE_SEGMENTS, DEFAULT_WIDTH, DEFAULT_HEIGHT, UPDATE_FIG_HEIGHT_FOR_MARGINS, GENE_NAME_cex, PLOT_PIPELINE_TRANSITIONS, DIFFERENTIATE_SAMPLE_GROUPS, COLOR_DIFFERENTIAL_TARGETS, MARK_SAMPLES, MARK_SAMPLES_COLORS, MARK_SAMPLES_ONLY_IN_CALLS, EXCLUDE_SAMPLES, MARK_INTERVALS, xlim=c(start, stop))
}

plot_XHMM_genes <- function(BASE_OUTPUT_NAME, xhmm_data, allTargsToGenes, binarySampleFeatures, SQ_THRESH, genes, PLOT_ONLY_PNG=TRUE, PLOT_LINE_SEGMENTS=TRUE, DEFAULT_WIDTH=14, DEFAULT_HEIGHT=7, UPDATE_FIG_HEIGHT_FOR_MARGINS=TRUE, GENE_NAME_cex=0.8, EXACT_GENE_MATCH=FALSE, PLOT_PIPELINE_TRANSITIONS=TRUE, DIFFERENTIATE_SAMPLE_GROUPS=NULL, COLOR_DIFFERENTIAL_TARGETS=TRUE, MARK_SAMPLES=NULL, MARK_SAMPLES_COLORS=DEFAULT_MARK_SAMPLES_COLORS, MARK_SAMPLES_ONLY_IN_CALLS=FALSE, EXCLUDE_SAMPLES=NULL, MARK_INTERVALS=NULL) {
	OUTPUT_NAME = paste(BASE_OUTPUT_NAME, ".genes_", paste(genes, collapse="_"), sep="")

	genesArray = genes
	if (EXACT_GENE_MATCH) {
		genesArray = paste("^", genes, "$", sep="")
	}
	genesGrepExpr = paste(genesArray, collapse="|")
	
	foundTargsToGenes = allTargsToGenes[grep(genesGrepExpr, allTargsToGenes, perl=TRUE)]
	foundTargs = names(foundTargsToGenes)
	
	removedTargs = setdiff(foundTargs, colnames(xhmm_data[["PCA_NORM_Z_SCORES"]]))
	if (length(removedTargs) > 0) {
		writeLines(paste(length(removedTargs), " ", genes, " targets were filtered out: ", paste(removedTargs, collapse=", "), sep=""))
	}
	
	foundTargs = intersect(foundTargs, colnames(xhmm_data[["PCA_NORM_Z_SCORES"]]))
	if (length(foundTargs) == 0) {
		writeLines(paste("No targets fully processed for genes: ", paste(genes, collapse=", "), sep=""))
		return(1)
	}
	foundTargsToGenes = foundTargsToGenes[foundTargs]

	foundTargs_chrBp1Bp2 = targetsToChrBp1Bp2(foundTargs)
	genesChr = unique(foundTargs_chrBp1Bp2[["chr"]])
	if (length(genesChr) != 1) {
		stop(paste("Genes: ", paste(genes, collapse=", "), " does not lie on a SINGLE chromosome [but rather: ", paste(genesChr, collapse=", "), "]", sep=""))
	}

	plot_XHMM_region(OUTPUT_NAME, xhmm_data, allTargsToGenes, binarySampleFeatures, SQ_THRESH, genesChr[1], min(foundTargs_chrBp1Bp2[["bp1"]]), max(foundTargs_chrBp1Bp2[["bp2"]]), PLOT_ONLY_PNG, PLOT_LINE_SEGMENTS, DEFAULT_WIDTH, DEFAULT_HEIGHT, UPDATE_FIG_HEIGHT_FOR_MARGINS, GENE_NAME_cex, PLOT_PIPELINE_TRANSITIONS, DIFFERENTIATE_SAMPLE_GROUPS, COLOR_DIFFERENTIAL_TARGETS, MARK_SAMPLES, MARK_SAMPLES_COLORS, MARK_SAMPLES_ONLY_IN_CALLS, EXCLUDE_SAMPLES, MARK_INTERVALS, APPEND_REGION_NAME=FALSE)
}

plot_XHMM_targets <- function(BASE_OUTPUT_NAME, xhmm_data, allTargsToGenes, binarySampleFeatures, SQ_THRESH, startTarg, stopTarg, PLOT_ONLY_PNG=TRUE, PLOT_LINE_SEGMENTS=TRUE, DEFAULT_WIDTH=14, DEFAULT_HEIGHT=7, UPDATE_FIG_HEIGHT_FOR_MARGINS=TRUE, GENE_NAME_cex=0.8, PLOT_PIPELINE_TRANSITIONS=TRUE, DIFFERENTIATE_SAMPLE_GROUPS=NULL, COLOR_DIFFERENTIAL_TARGETS=TRUE, MARK_SAMPLES=NULL, MARK_SAMPLES_COLORS=DEFAULT_MARK_SAMPLES_COLORS, MARK_SAMPLES_ONLY_IN_CALLS=FALSE, EXCLUDE_SAMPLES=NULL, MARK_INTERVALS=NULL, APPEND_REGION_NAME=TRUE) {
	chrBp1Bp2 = targetsToChrBp1Bp2(colnames(xhmm_data[["PCA_NORM_Z_SCORES"]]))
	chr = chrBp1Bp2[["chr"]]

	if (startTarg < 1 || stopTarg > length(chr) || startTarg > stopTarg || chr[startTarg] != chr[stopTarg]) {
		stop(paste("Invalid target start and stop: ", startTarg, "-", stopTarg, sep=""))
	}

	plot_XHMM_region(BASE_OUTPUT_NAME, xhmm_data, allTargsToGenes, binarySampleFeatures, SQ_THRESH, chr[startTarg], chrBp1Bp2[["bp1"]][startTarg], chrBp1Bp2[["bp2"]][stopTarg], PLOT_ONLY_PNG, PLOT_LINE_SEGMENTS, DEFAULT_WIDTH, DEFAULT_HEIGHT, UPDATE_FIG_HEIGHT_FOR_MARGINS, GENE_NAME_cex, PLOT_PIPELINE_TRANSITIONS, DIFFERENTIATE_SAMPLE_GROUPS, COLOR_DIFFERENTIAL_TARGETS, MARK_SAMPLES, MARK_SAMPLES_COLORS, MARK_SAMPLES_ONLY_IN_CALLS, EXCLUDE_SAMPLES, MARK_INTERVALS, APPEND_REGION_NAME)
}



#####################################################
# The implementation:
#####################################################

parseRegion <- function(chrom, start=NULL, stop=NULL) {
	if (is.null(start) || is.null(stop)) {
		parse1 = strsplit(chrom, ":")
		chrom = sapply(parse1, "[[", 1)

		parse2 = sapply(parse1, "[[", 2)
		bp1_bp2 = strsplit(parse2, "-")

		start = as.numeric(sapply(bp1_bp2, "[[", 1))
		stop = as.numeric(sapply(bp1_bp2, "[[", 2))
	}

	return(list(chrom=chrom, start=start, stop=stop))
}

regionToFinalInds <- function(xhmm_data, chrom, start=NULL, stop=NULL) {
	if (is.null(start) || is.null(stop)) {
		chrom_start_stop = parseRegion(chrom, start, stop)
		chrom = as.character(chrom_start_stop["chrom"])
		start = as.numeric(chrom_start_stop["start"])
		stop  = as.numeric(chrom_start_stop["stop"])
	}

	chrBp1Bp2 = targetsToChrBp1Bp2(colnames(xhmm_data[["PCA_NORM_Z_SCORES"]]))
	return(which(chrBp1Bp2[["chr"]] == chrom & chrBp1Bp2[["bp1"]] <= stop & start <= chrBp1Bp2[["bp2"]]))
}

regionToStartStopInds <- function(xhmm_data, chrom, start=NULL, stop=NULL, NUM_ADD_TARGS=2) {
	targFinalInds = regionToFinalInds(xhmm_data, chrom, start, stop)

	startTarg = min(targFinalInds)
	stopTarg = max(targFinalInds)

	chr = targetsToChrBp1Bp2(colnames(xhmm_data[["PCA_NORM_Z_SCORES"]]))[["chr"]]

	if (NUM_ADD_TARGS > 0) {
		potentialStartTargs = max(1, startTarg - NUM_ADD_TARGS):startTarg
		startTarg = min(potentialStartTargs[chr[potentialStartTargs] == chr[startTarg]])

		potentialStopTargs = stopTarg:min(length(chr), stopTarg + NUM_ADD_TARGS)
		stopTarg = max(potentialStopTargs[chr[potentialStopTargs] == chr[stopTarg]])
	}

	return(c(startTarg=startTarg, stopTarg=stopTarg))
}


DELETION = 1
DIPLOID = 2
DUPLICATION = 3

CNVnameToIndex <- function(CNV_names) {
	CNV_indices = CNV_names
	CNV_indices = 0

	CNV_indices[CNV_names == "DEL"] = DELETION
	CNV_indices[CNV_names == "DIP"] = DIPLOID
	CNV_indices[CNV_names == "DUP"] = DUPLICATION

	return(CNV_indices)
}

col2rgbString <- function(color, alphaFactor=1) {
	MAX_VAL = 255
	return(apply(col2rgb(color), 2, function(x) rgb(x[1], x[2], x[3], alpha=alphaFactor*MAX_VAL, maxColorValue=MAX_VAL)))
}

CNVtoColor <- function(dataVals) {
	colors = dataVals
	colors[] = col2rgbString("gray")

	colors[abs(dataVals) < DIPLOID] = col2rgbString("red")
	colors[abs(dataVals) > DIPLOID] = col2rgbString("green")

	# Below-threshold calls:
	belowThreshInds = which(sign(dataVals) == -1)
	colors[belowThreshInds] = col2rgbString(colors[belowThreshInds], 0.2)

	return(colors)
}

CNVtoWidth <- function(dataVals) {
	widths = dataVals
	widths[] = 2

	widths[abs(dataVals) != DIPLOID] = 4

	return(widths)
}

startNewRow <- function(intervInd, interv) {
	return(intervInd == 1 || length(which(names(interv) == "SAME_ROW")) == 0 || interv[["SAME_ROW"]] != 1)
}

plot_XHMM_targets_helper <- function(BASE_OUTPUT_NAME, xhmm_data, allTargsToGenes, binarySampleFeatures, SQ_THRESH, targs, PLOT_ONLY_PNG=TRUE, PLOT_LINE_SEGMENTS=TRUE, DEFAULT_WIDTH=14, DEFAULT_HEIGHT=7, UPDATE_FIG_HEIGHT_FOR_MARGINS=TRUE, GENE_NAME_cex=0.8, PLOT_PIPELINE_TRANSITIONS=TRUE, DIFFERENTIATE_SAMPLE_GROUPS=NULL, COLOR_DIFFERENTIAL_TARGETS=TRUE, MARK_SAMPLES=NULL, MARK_SAMPLES_COLORS=DEFAULT_MARK_SAMPLES_COLORS, MARK_SAMPLES_ONLY_IN_CALLS=FALSE, EXCLUDE_SAMPLES=NULL, MARK_INTERVALS=NULL, xlim=NULL) {
	graphics.off()

	ALL_SAMPLES = rownames(xhmm_data[["PCA_NORM_Z_SCORES"]])
	if (!is.null(EXCLUDE_SAMPLES)) {
		ALL_SAMPLES = setdiff(ALL_SAMPLES, EXCLUDE_SAMPLES)
	}

	# Check which targets are present in the final matrix:
	ALL_TARGS_TO_FINAL_INDS = 1:ncol(xhmm_data[["PCA_NORM_Z_SCORES"]])
	names(ALL_TARGS_TO_FINAL_INDS) = colnames(xhmm_data[["PCA_NORM_Z_SCORES"]])
	targFinalInds = as.numeric(ALL_TARGS_TO_FINAL_INDS[targs])

	notNAinds = which(!is.na(targFinalInds))
	targFinalInds = targFinalInds[notNAinds]
	targs = targs[notNAinds]

	targs_chrBp1Bp2 = targetsToChrBp1Bp2(targs)

	if (is.null(allTargsToGenes)) {
		allTargsToGenes = list()
	}
	targsToGenes = allTargsToGenes[targs]
	targs_genes = setdiff(unique(as.character(unlist(targsToGenes))), ".")

	MAX_GENES_TITLE = 6
	targsAllGenesString = paste(targs_genes[1:min(MAX_GENES_TITLE, length(targs_genes))], collapse=";")
	targsToGeneString = sapply(targsToGenes, function(x) paste(x, collapse=";"))

	if (length(targs_genes) > 0) {
		TMP_prepend_MARK_INTERVALS = sapply(targs_genes, function(singleGene) {gene_inds = which(as.logical(lapply(targsToGenes, function(x) any(x == singleGene)))); gene_chrBp1Bp2 = targetsToChrBp1Bp2(targs[gene_inds]); return(list(min(gene_chrBp1Bp2[["bp1"]]), max(gene_chrBp1Bp2[["bp2"]]), cex=GENE_NAME_cex, col="dodgerblue2", ROW_IND=max(as.numeric(sapply(targsToGenes[gene_inds], function(x) which(sort(x) == singleGene))), na.rm=TRUE)))}, simplify=FALSE, USE.NAMES=TRUE)

		NUM_GENE_TRACKS = max(sapply(targsToGenes, function(x) length(x[x != "."])))
		GENE_TRACKS = list()
		for (i in 1:NUM_GENE_TRACKS) {GENE_TRACKS[[i]]=list()}
	
		for (i in 1:length(TMP_prepend_MARK_INTERVALS)) {
			interv = TMP_prepend_MARK_INTERVALS[[i]]
			geneName = names(TMP_prepend_MARK_INTERVALS)[i]
			trackInd = interv[["ROW_IND"]]
			GENE_TRACKS[[trackInd]][[geneName]] = interv
		}
	
		# Put the gene tracks in order and mark where they should start:
		prepend_MARK_INTERVALS = list()
		for (i in 1:NUM_GENE_TRACKS) {
			iTracks = GENE_TRACKS[[i]]
			for (j in 1:length(iTracks)) {
				geneName = names(iTracks[j])
				geneTrack = iTracks[[j]]
				if (j > 1) {
					geneTrack[["SAME_ROW"]] = 1
				}
				prepend_MARK_INTERVALS[[geneName]] = geneTrack
			}
		}
	
		if (is.null(MARK_INTERVALS)) {
			MARK_INTERVALS = list()
		}
		MARK_INTERVALS = c(prepend_MARK_INTERVALS, MARK_INTERVALS)
	}

	ALL_CNV = xhmm_data[["XCNV_CALLS"]]

	OVERLAP_TARG_CNV_INDS = which(sapply(strsplit(ALL_CNV[, "TARGETS"], "\\.\\."), function(x) {x=as.numeric(x); length(intersect(x[1]:x[2], targFinalInds)) > 0}))
	OVERLAP_TARG_CNV = ALL_CNV[OVERLAP_TARG_CNV_INDS, ]
	OVERLAP_TARG_CNV = OVERLAP_TARG_CNV[ OVERLAP_TARG_CNV[, "SAMPLE"] %in% ALL_SAMPLES , ]

	HAVE_CNV_SAMPLES = intersect(unique(OVERLAP_TARG_CNV[, "SAMPLE"]), ALL_SAMPLES)
	HAVE_ABOVE_THRESH_CNV_SAMPLES = intersect(unique(OVERLAP_TARG_CNV[OVERLAP_TARG_CNV[, "Q_SOME"] >= SQ_THRESH, "SAMPLE"]), HAVE_CNV_SAMPLES)
	HAVE_CNV_SAMPLES = c(setdiff(HAVE_CNV_SAMPLES, HAVE_ABOVE_THRESH_CNV_SAMPLES), HAVE_ABOVE_THRESH_CNV_SAMPLES)

	CNV_calls = matrix(DIPLOID, nrow=length(ALL_SAMPLES), ncol=length(targs), dimnames=list(ALL_SAMPLES, targs))
	sample_extent_CNV_SQs = apply(OVERLAP_TARG_CNV, 1, function(x) {t1_t2=strsplit(x["TARGETS"], "\\.\\."); list(as.character(x["SAMPLE"]), as.numeric(unlist(t1_t2)), CNVnameToIndex(as.character(x["CNV"])), as.numeric(x["Q_SOME"]))})
	for (i in 1:length(sample_extent_CNV_SQs)) {
		sample_extent_CNV_SQ = sample_extent_CNV_SQs[[i]]
		sample = sample_extent_CNV_SQ[[1]]
		extent = sample_extent_CNV_SQ[[2]]
		CNV = sample_extent_CNV_SQ[[3]]
		SQ = sample_extent_CNV_SQ[[4]]

		multFactor = 1
		if (SQ < SQ_THRESH) {
			multFactor = -1
		}

		local_CNV_inds = which(targFinalInds %in% extent[1]:extent[2])
		CNV_calls[sample, local_CNV_inds] = multFactor * CNV
	}

	if (0) {
		haveCNVtargLocalInds = which(apply(abs(CNV_calls), 2, function(x) any(x != DIPLOID)))
		haveCNVtargs = targs[haveCNVtargLocalInds]
		writeLines(paste("Targets with some CNV: ", paste(haveCNVtargs, collapse=", "), sep=""))
		interestingStart = min(targs_chrBp1Bp2[["bp1"]][haveCNVtargLocalInds])
		interestingStop  = max(targs_chrBp1Bp2[["bp2"]][haveCNVtargLocalInds])
		writeLines(paste("Sub-range: ", targs_chrBp1Bp2[["chr"]][haveCNVtargLocalInds[1]], ":", interestingStart, "-", interestingStop, sep=""))
		return(-1)
	}
	
	startBP = min(targs_chrBp1Bp2[["bp1"]])
	stopBP = max(targs_chrBp1Bp2[["bp2"]])

	if (!is.null(xlim) && length(xlim) == 2) {
		startBP = min(startBP, xlim[1])
		stopBP  = max(stopBP,  xlim[2])
	}
	extent = startBP:stopBP

	logarithm = round(log10(length(extent)), 0)
	title1 = paste(targsAllGenesString, " [", length(targFinalInds), " targets in ", sep="")

	addChrString = ""
	if (!xhmm_data[["CHROMOSOMES_START_WITH_CHR"]]) {
		addChrString = "chr"
	}
	title2 = paste(" bases, ", addChrString, unique(targs_chrBp1Bp2[["chr"]]), ":", format(startBP, scientific=FALSE), "-", format(stopBP, scientific=FALSE), "]", sep="")

	title = eval(substitute(expression(paste(title1, 10 ^ logarithm, title2, sep="")), list(logarithm=logarithm, title1=title1, title2=title2)))
	
	PLOT_LAST = HAVE_CNV_SAMPLES
	if (!is.null(MARK_SAMPLES)) {
		HAVE_MARK_SAMPLES = intersect(MARK_SAMPLES, ALL_SAMPLES)
		PLOT_LAST = union(PLOT_LAST, HAVE_MARK_SAMPLES)
		PLOT_LAST = c(setdiff(PLOT_LAST, HAVE_MARK_SAMPLES), HAVE_MARK_SAMPLES)

		if (length(MARK_SAMPLES) <= length(MARK_SAMPLES_COLORS)) {
			MARK_SAMPLES_COLORS = MARK_SAMPLES_COLORS[1:length(MARK_SAMPLES)]
		}
		else {
			MARK_SAMPLES_COLORS = heat.colors(length(MARK_SAMPLES))
		}

		names(MARK_SAMPLES_COLORS) = MARK_SAMPLES

		MARK_SAMPLES_pch = "o"
		MARK_SAMPLES_cex = 1.5
	}

	USE_SAMPLES = c(setdiff(ALL_SAMPLES, PLOT_LAST), PLOT_LAST)


	if (!is.null(DIFFERENTIATE_SAMPLE_GROUPS)) {
		CONTROL_SAMPLES = intersect(USE_SAMPLES, DIFFERENTIATE_SAMPLE_GROUPS[[1]])
		CASE_SAMPLES    = intersect(USE_SAMPLES, DIFFERENTIATE_SAMPLE_GROUPS[[2]])
	}

	P_VAL_THRESH = 0.05




	READ_DEPTH_PIPELINE_MAT_NAMES = c("RD", "filtered_centered_RD", "PCA_NORMALIZED", "PCA_NORM_Z_SCORES", "DEL_POST", "DUP_POST", "DIP_POST")
	READ_DEPTH_PIPELINE_SHORT_NAMES = c("read depth", "Centering", "PCA", "Sample z-score", "P(deletion)", "P(duplication)", "P(diploid)")
	READ_DEPTH_PIPELINE_FULL_NAMES = c("Unnormalized read depth", "Target-centered read depth", "PCA-normalized read depth", "Sample Z-score of PCA-normalized read depth", "P(deletion)", "P(duplication)", "P(diploid)")

	PROB_LIMS = c(0, 1)
	YLIMS = list(NULL, NULL, NULL, NULL, PROB_LIMS, PROB_LIMS, PROB_LIMS)


	READ_DEPTH_PIPELINE = listOfNulls(length(READ_DEPTH_PIPELINE_MAT_NAMES))
	for (ind in 1:length(READ_DEPTH_PIPELINE_MAT_NAMES)) {
		mat = xhmm_data[[ READ_DEPTH_PIPELINE_MAT_NAMES[ind] ]]
		if (!is.null(mat)) {
			READ_DEPTH_PIPELINE[[ind]] = as.matrix(mat[USE_SAMPLES, targs])
		}
	}

	NORMALIZED_RD_IND = 4
	PIPELINE_PAIRS = rbind(c(1, 2), c(2, 3), c(3, 4))

	DEFAULT_TEXT_CEX = 3

	HEIGHT = DEFAULT_HEIGHT

	# margins in the form: c(bottom, left, top, right):
	#margins = par("mar")
	margins = c(5, 4, 4, 2) + 0.1  # R's default

	if (!is.null(MARK_INTERVALS) && length(MARK_INTERVALS) > 0) {
		MARGIN_TO_ADD = 0
		maxCex = -1
		relMarginsToAdd = c()
		for (intervInd in 1:length(MARK_INTERVALS)) {
			interv = MARK_INTERVALS[[intervInd]]

			cex = DEFAULT_TEXT_CEX
			cexInds = which(names(interv) == "cex")
			if (length(cexInds) > 0) {
				cex = interv[[cexInds[1]]]
			}

			if (startNewRow(intervInd, interv)) {
				if (maxCex != -1) {
					marginToAdd = (0.15 + (maxCex * 1.2))
					MARGIN_TO_ADD = MARGIN_TO_ADD + marginToAdd
					relMarginsToAdd = c(relMarginsToAdd, marginToAdd)
				}
				maxCex = cex
			}

			maxCex = max(maxCex, cex)
		}
		if (maxCex != -1) {
			marginToAdd = (0.15 + (maxCex * 1.2))
			MARGIN_TO_ADD = MARGIN_TO_ADD + marginToAdd
			relMarginsToAdd = c(relMarginsToAdd, marginToAdd)
		}
		relMarginsToAdd = relMarginsToAdd / sum(relMarginsToAdd)

		TOTAL_MARGIN = margins[1] + MARGIN_TO_ADD
		margins[1] = TOTAL_MARGIN
		TOTAL_MARGIN_TO_USE = 0.8 * TOTAL_MARGIN

		FRAC_MARG_FOR_TICKS = 0.12

		if (UPDATE_FIG_HEIGHT_FOR_MARGINS) {
			HEIGHT = DEFAULT_HEIGHT * (1 + MARGIN_TO_ADD/TOTAL_MARGIN)
		}
	}

	if (!PLOT_ONLY_PNG) {
		PIPELINE_INDS = 1:length(READ_DEPTH_PIPELINE)
		OUTPUT_FILE = paste(BASE_OUTPUT_NAME, ".pdf", sep="")
		pdf(OUTPUT_FILE, width=DEFAULT_WIDTH, height=HEIGHT)
	}
	else {
		PIPELINE_INDS = NORMALIZED_RD_IND
		OUTPUT_FILE = paste(BASE_OUTPUT_NAME, ".png", sep="")
		png(OUTPUT_FILE, width=DEFAULT_WIDTH, height=HEIGHT, units="in", res=150)
	}
	writeLines(paste("Output: ", OUTPUT_FILE, sep=""))

	# TO PERMIT PLOTTING "INTERVALS" LINE SEGMENTS BELOW THE PLOT:
	par(xpd=TRUE)
	par(mar=margins)

	for (typeInd in PIPELINE_INDS) {
		data = READ_DEPTH_PIPELINE[[typeInd]]
		dataName = READ_DEPTH_PIPELINE_FULL_NAMES[typeInd]
		
		ylim = YLIMS[[typeInd]]
		if (is.null(ylim)) {
			ylim = c(min(data, na.rm=TRUE), max(data, na.rm=TRUE))
		}

		xlab = "Genome locus"
		if (!is.null(MARK_INTERVALS)) {
			xlab=""
		}

		oldScipen = getOption("scipen")
		options(scipen = 50)
		plot(extent, rep(0, length(extent)), type="n", ylim=ylim, xlab=xlab, ylab=dataName, main=title)
		options(scipen = oldScipen)
		
		# Plot the target intervals:
		yLimScale = ylim[2] - ylim[1]
		plotYval = ylim[1] - 0.025 * yLimScale
		for (i in 1:length(targs)) {
			x = c(targs_chrBp1Bp2[["bp1"]][i], targs_chrBp1Bp2[["bp2"]][i])
			y = rep(plotYval, 2)

			type = "l"
			lty = "solid"
			pch = "."
			col = "black"

			if (!is.null(DIFFERENTIATE_SAMPLE_GROUPS)) {
				controlData = data[CONTROL_SAMPLES, i]
				caseData = data[CASE_SAMPLES, i]
	
				if (length(which(!is.na(controlData))) > 0 && length(which(!is.na(caseData))) > 0) {
					caseControlTtest = t.test(controlData, caseData)
					if (caseControlTtest$p.val <= P_VAL_THRESH) {
						if (0) {
							writeLines(paste("SIGNIFICANT groups difference for target ", targs[i], " in ", dataName, sep=""))
							print(caseControlTtest)
						}
		
						type = "o"
						lty = "solid"
		
						if (COLOR_DIFFERENTIAL_TARGETS) {
							controlMean = mean(data[CONTROL_SAMPLES, i])
							caseMean = mean(data[CASE_SAMPLES, i])
							if (caseMean < controlMean) {
								pch = "o"
								col = "darkgoldenrod2"
							}
							else {
								pch = "+"
								col = "blue"
							}
						}
					}
				}
			}

			points(x, y, type=type, col=col, lwd=5, pch=pch, lty=lty)
		}

		if (!is.null(MARK_INTERVALS) && length(MARK_INTERVALS) > 0) {
			addIndex = 1
			for (intervInd in 1:length(MARK_INTERVALS)) {
				intervName = names(MARK_INTERVALS)[intervInd]
				interv = MARK_INTERVALS[[intervInd]]

				intervStart = as.numeric(interv[[1]])
				intervStop = as.numeric(interv[[2]])
				intervCenter = (intervStart + intervStop) / 2

				col = "black"
				colInds = which(names(interv) == "col")
				if (length(colInds) > 0) {
					col = interv[[colInds[1]]]
				}

				cex = DEFAULT_TEXT_CEX
				cexInds = which(names(interv) == "cex")
				if (length(cexInds) > 0) {
					cex = interv[[cexInds[1]]]
				}

				if (startNewRow(intervInd, interv)) {
					sumRelMargin = sum(relMarginsToAdd[1:addIndex])
					linePos = sumRelMargin * TOTAL_MARGIN_TO_USE
					fracMargUp = 1 - (linePos / TOTAL_MARGIN) - FRAC_MARG_FOR_TICKS

					addIndex = addIndex + 1
				}

				# Convert the plot's 0 to normalized inner region coordinates, 
				# go up by the relative fraction of this "track", 
				# and then convert to user coordinates (as required by segments):
				plotYval = grconvertY(grconvertY(0, from = "npc", to = "nic") * fracMargUp, from = "nic", to = "user")
				segments(intervStart, plotYval, intervStop, plotYval, col=col, lwd=8)

				text(intervCenter, plotYval, intervName, col=col, adj=c(0.5, -0.75), cex=cex)
				#mtext(intervName, side=1, line=linePos, at=intervCenter, col=col, cex=cex, padj=0)

				lineTextsInd = which(names(interv) == "lineTexts")
				for (ind in lineTextsInd) {
					textPositions = interv[[ind]]
					for (tInd in 1:length(textPositions)) {
						text(textPositions[tInd], plotYval, names(textPositions)[tInd], col="black", adj=c(0.5, 0.45), cex=1.8)
					}
				}
			}
		}
		
		targetsPos = as.vector(matrix(-1, nrow=2*length(targs), ncol=1))
		STOP_TARGS = 2*(1:length(targs))
		START_TARGS = STOP_TARGS - 1
		targetsPos[START_TARGS] = targs_chrBp1Bp2[["bp1"]]
		targetsPos[STOP_TARGS] = targs_chrBp1Bp2[["bp2"]]

		middleTargetsPos = apply(rbind(targetsPos[START_TARGS], targetsPos[STOP_TARGS]), 2, mean)
		
		for (s in 1:length(USE_SAMPLES)) {
			samp = USE_SAMPLES[s]
		
			y = as.vector(matrix(0, nrow=2*length(targs), ncol=1))
			y[START_TARGS] = data[s, ]
			y[STOP_TARGS] =  data[s, ]
		
			singleCallVec = as.vector(CNV_calls[samp, targs])
			callVec = as.vector(matrix(0, nrow=2*length(targs), ncol=1))
			callVec[START_TARGS] = singleCallVec
			callVec[STOP_TARGS] =  singleCallVec
	
	
			# Find indices where this call and the subsequent one are the same:
			consecCallVec = callVec
			consecCallVec[STOP_TARGS] = DIPLOID

			# All indices BUT the ends of intervals are, by definition, in the "interior" of an interval:
			startEndIntervals = calcSegmentIntervals(list(singleCallVec))
			endInds = startEndIntervals[, 2]
			revertInds = setdiff(1:length(singleCallVec), endInds)

			consecCallVec[STOP_TARGS[revertInds]] = singleCallVec[revertInds]
		
			lineColVec = CNVtoColor(consecCallVec)
	
			widthVec = CNVtoWidth(consecCallVec)
			widthVec[STOP_TARGS] = widthVec[STOP_TARGS] / 3
	
			ltyVec = as.vector(matrix(0, nrow=2*length(targs), ncol=1))
			ltyVec[] = "solid"  # Used by default, and for "cases" (when there are two groups)
			if (!is.null(DIFFERENTIATE_SAMPLE_GROUPS)) {
				if (samp %in% CONTROL_SAMPLES) {
					ltyVec[] = "dashed"
				}
				else if (!(samp %in% CASE_SAMPLES)) {
					# Not in "cases" or "controls", so mark as having no group:
					ltyVec[] = "twodash"
				}
			}

			if (!PLOT_LINE_SEGMENTS) {
				stopIndsInCnv = which(consecCallVec != DIPLOID)
				ltyVec[setdiff(STOP_TARGS, stopIndsInCnv)] = "blank"
			}

			color.scale.lines(targetsPos, y, col=lineColVec, lwd=widthVec, lty=ltyVec)

			if (!is.null(MARK_SAMPLES) && samp %in% MARK_SAMPLES) {
				markInds = 1:length(targs)
				if (MARK_SAMPLES_ONLY_IN_CALLS) {
					markInds = (singleCallVec != DIPLOID)
				}
				points(middleTargetsPos[markInds], data[s, markInds], type="p", col=MARK_SAMPLES_COLORS[samp], pch=MARK_SAMPLES_pch, cex=MARK_SAMPLES_cex)
			}
		}

		if (!is.null(MARK_SAMPLES)) {
			REPLACE_PREFIX = "::REPLACE::"

			MARK_SAMPLES_LEGEND = MARK_SAMPLES
			if (!is.null(names(MARK_SAMPLES))) {
				MARK_SAMPLES_LEGEND = as.character(sapply(1:length(MARK_SAMPLES), function(i) {str=MARK_SAMPLES[i]; nm=names(MARK_SAMPLES)[i]; if (nm != "") {matches=grep(REPLACE_PREFIX, nm, value=FALSE); if (length(matches) == 1 && matches == 1) {str = substring(nm, nchar(REPLACE_PREFIX)+1)} else {str = paste(str, " (", nm, ")", sep="")}}; return(str)}))
			}

			legend("topright", MARK_SAMPLES_LEGEND, inset=c(-0.02, -0.05), xpd=NA, col=MARK_SAMPLES_COLORS, pch=MARK_SAMPLES_pch, pt.cex=MARK_SAMPLES_cex, bg="white")
		}
	}
	dev.off()

	if (PLOT_ONLY_PNG || !PLOT_PIPELINE_TRANSITIONS) {
		return(1)
	}
	




	# Plot the pipeline's processing of each target with respect to the various sample labels:
	if (!is.null(binarySampleFeatures)) {
		LABEL = binarySampleFeatures
	}
	else {
		LABEL = list(All=array(0, dim=length(ALL_SAMPLES), dimnames=list(ALL_SAMPLES)))
	}

	LABEL_NAME = names(LABEL)
	SHOW_LABELS = replicate(length(LABEL), NULL)
	
	
	pch = "."
	cex = 4

	NUM_PLOTS = length(LABEL) * (nrow(PIPELINE_PAIRS) + 1)
	LAYOUT_MATRIX = matrix(1:NUM_PLOTS, nrow=(nrow(PIPELINE_PAIRS) + 1), byrow=FALSE)
	LAYOUT_MATRIX = rbind(1, LAYOUT_MATRIX + 1)
	
	LAYOUT_HEIGHTS = rep(1, nrow(LAYOUT_MATRIX))

	# title row:
	LAYOUT_HEIGHTS[1] = 0.1

	# legend row:
	LAYOUT_HEIGHTS[2] = 0.07 * max(sapply(LABEL, function(x) length(unique(as.numeric(x)))))
	

	PIPELINE_OUTPUT_FILE = paste(BASE_OUTPUT_NAME, ".pipeline.pdf", sep="")
	pdf(PIPELINE_OUTPUT_FILE, width=10, height=25)
	writeLines(paste("Output: ", PIPELINE_OUTPUT_FILE, sep=""))
	
	#
	if (0) {
		nf = layout(LAYOUT_MATRIX, heights = LAYOUT_HEIGHTS)
		layout.show(nf)
	}
	#
	
	for (k in 1:length(targs)) {
		layout(LAYOUT_MATRIX, heights = LAYOUT_HEIGHTS)
	
		# margins in the form: c(bottom, left, top, right):
		par(mar=c(0, 0, 2.2, 0), xpd=TRUE)
		plot.new()
		title = paste("Target ", targs[k], " [", targsToGeneString[k], "]", sep="")
		title(main=list(title, cex=2.5))
	
		for (j in 1:length(LABEL)) {
			labels = LABEL[[j]]
			unique_labels = sort(unique(labels))
			labels_to_numbers = array(1:length(unique_labels), dim=length(unique_labels), dimnames=list(unique_labels))
			sampleColors = array(-1, dim=length(labels), dimnames=list(names(labels)))
			for (l in 1:length(labels_to_numbers)) {sampleColors[labels == unique_labels[l]] = l}

			legLabels = SHOW_LABELS[[j]]
			if (!is.null(legLabels)) {
				legLabels = legLabels[unique_labels]
			}
			else {
				legLabels = as.character(unique_labels)
			}
	
			# margins in the form: c(bottom, left, top, right):
			par(mar=c(2, 1, 2, 1), xpd=TRUE)
			plot.new()
			title(LABEL_NAME[[j]])
			legend("center", legLabels, pch=pch, col=as.numeric(labels_to_numbers), pt.cex=cex, xjust=0.5, yjust=0.5)
	
			for (i in 1:nrow(PIPELINE_PAIRS)) {
				i1 = PIPELINE_PAIRS[i, 1]
				i2 = PIPELINE_PAIRS[i, 2]

				plotSamples = intersect(names(sampleColors), intersect(rownames(READ_DEPTH_PIPELINE[[i1]]), rownames(READ_DEPTH_PIPELINE[[i2]])))
	
				# margins in the form: c(bottom, left, top, right):
				par(mar=c(4.1, 4.1, 1, 0.5), xpd=TRUE)
				plot(as.numeric(READ_DEPTH_PIPELINE[[i1]][plotSamples, k]), as.numeric(READ_DEPTH_PIPELINE[[i2]][plotSamples, k]), col=sampleColors[plotSamples], pch=pch, xlab=READ_DEPTH_PIPELINE_SHORT_NAMES[i1], ylab=READ_DEPTH_PIPELINE_SHORT_NAMES[i2], cex=cex)
			}
		}
	}
	
	dev.off()
}
