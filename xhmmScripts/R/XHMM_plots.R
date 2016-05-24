XHMM_plots <- function(PLOT_PATH, JOB_PREFICES, JOB_TARGETS_TO_GENES=NULL, SAMPLE_FEATURES=NULL, SQ_THRESH=60, NUM_ADD_TARGS=2, PLOT_READ_DEPTHS=TRUE, PLOT_PC_CORRS=TRUE, PLOT_ALL_CNVS=TRUE, USE_XCNV_TO_PLOT=NULL, INCLUDE_PEDIGREE_SAMPLES=NULL, PLOT_ONLY_PNG=TRUE, LIMIT_MEMORY=FALSE) {
	JOBS = names(JOB_PREFICES)
	if (is.null(JOBS)) {
		JOBS = paste("dataset_", as.character(1:length(JOB_PREFICES)), sep="")
		names(JOB_PREFICES) = JOBS

		if (!is.null(JOB_TARGETS_TO_GENES)) {
			if (length(JOB_TARGETS_TO_GENES) != length(JOBS) || !is.null(names(JOB_TARGETS_TO_GENES))) {
				stop("If JOB_PREFICES is not named, then JOB_TARGETS_TO_GENES must also be un-named and of the same length")
			}
			names(JOB_TARGETS_TO_GENES) = JOBS
		}
	}
	
	SAMPLE_FEATURE_NAMES = names(SAMPLE_FEATURES)

	READ_DEPTH_PIPELINE = c("RD", "filtered_centered_RD", "PCA_NORMALIZED", "PCA_NORM_Z_SCORES")
	READ_DEPTH_PIPELINE_SHORT_NAMES = c("DoC read-depth", "Target centering", "PCA", "Sample z-score")
	READ_DEPTH_PIPELINE_FULL_NAMES = c("Unnormalized read depth", "Target-centered read depth", "PCA-normalized read depth", "Sample Z-score of PCA-normalized read depth")

	EXCLUDE_LARGE_MATRICES = c()
	if (LIMIT_MEMORY) {
		EXCLUDE_LARGE_MATRICES = c("filtered_centered_RD", "PCA_NORMALIZED", "PCA_NORM_Z_SCORES", "RD_SAME_FILTERED", "DIP_POST", "DEL_POST", "DUP_POST")
		if (PLOT_READ_DEPTHS) {
			# Keep in for plotting standard deviation:
			EXCLUDE_LARGE_MATRICES = setdiff(EXCLUDE_LARGE_MATRICES, c("PCA_NORMALIZED"))
		}

		if (PLOT_ALL_CNVS) {
			# Keep in for plotting CNV:
			EXCLUDE_LARGE_MATRICES = setdiff(EXCLUDE_LARGE_MATRICES, c("PCA_NORM_Z_SCORES"))
			if (!PLOT_ONLY_PNG) {
				EXCLUDE_LARGE_MATRICES = setdiff(EXCLUDE_LARGE_MATRICES, c("filtered_centered_RD", "PCA_NORMALIZED", "DIP_POST", "DEL_POST", "DUP_POST"))
			}
		}

		names(READ_DEPTH_PIPELINE_SHORT_NAMES) = READ_DEPTH_PIPELINE
		names(READ_DEPTH_PIPELINE_FULL_NAMES) = READ_DEPTH_PIPELINE

		READ_DEPTH_PIPELINE = setdiff(READ_DEPTH_PIPELINE, EXCLUDE_LARGE_MATRICES)
		READ_DEPTH_PIPELINE_SHORT_NAMES = READ_DEPTH_PIPELINE_SHORT_NAMES[READ_DEPTH_PIPELINE]
		READ_DEPTH_PIPELINE_FULL_NAMES = READ_DEPTH_PIPELINE_FULL_NAMES[READ_DEPTH_PIPELINE]
	}
	
	TARGET_FEATURE_NAMES = c("GC", "Repeat-masked", "Target size")
	
	
	xhmm_data = list()
	for (job in JOBS) {
		xhmm_data[[job]] = loadXHMMdata(JOB_PREFICES[job], EXCLUDE_LARGE_MATRICES=EXCLUDE_LARGE_MATRICES)
	}

	# Load genes for each target:
	if (!is.null(JOB_TARGETS_TO_GENES)) {
		for (job in JOBS) {
			if (!is.null(JOB_TARGETS_TO_GENES[job])) {
				writeLines(paste("Loading genes for ", job, " targets", sep=""))
				xhmm_data[[job]][["TARGETS_TO_GENES"]] = loadTargetsToGenes(JOB_TARGETS_TO_GENES[job], CHROMOSOMES_START_WITH_CHR=xhmm_data[[job]][["CHROMOSOMES_START_WITH_CHR"]])
			}
		}
	}

	# Copy sample data into per-job lists:
	for (job in JOBS) {
		for (sampVal in SAMPLE_FEATURE_NAMES) {
			xhmm_data[[job]][[sampVal]] = SAMPLE_FEATURES[[sampVal]]
		}
	}
	
	
	mean.target = list()
	mean.sample = list()
	for (rd in READ_DEPTH_PIPELINE) {
		for (job in JOBS) {
			mean.sample[[job]][[rd]] = rowMeans(xhmm_data[[job]][[rd]])
			mean.target[[job]][[rd]] = colMeans(xhmm_data[[job]][[rd]])
		}
	}
	
	if (PLOT_READ_DEPTHS) {
		pdf(paste(PLOT_PATH, "/mean_target_coverage.pdf", sep=""))
		for (job in JOBS) {
			for (i in 1:length(READ_DEPTH_PIPELINE)) {
				rd = READ_DEPTH_PIPELINE[i]
				rdName = READ_DEPTH_PIPELINE_SHORT_NAMES[i]
				hist(mean.target[[job]][[rd]], breaks=400, xlab="Mean target coverage", ylab="Frequency", main=paste(job, ": ", rdName, sep=""))
			}
		}
		dev.off()
		
		pdf(paste(PLOT_PATH, "/mean_sample_coverage.pdf", sep=""))
		for (job in JOBS) {
			for (i in 1:length(READ_DEPTH_PIPELINE)) {
				rd = READ_DEPTH_PIPELINE[i]
				rdName = READ_DEPTH_PIPELINE_SHORT_NAMES[i]
				hist(mean.sample[[job]][[rd]], breaks=400, xlab="Mean sample coverage", ylab="Frequency", main=paste(job, ": ", rdName, sep=""))
			}
		}
		dev.off()
	}
	
	
	
	# Add in per-target and per-sample mean read depths to correlate with PC and loadings:
	TARG_DEPTH = "Target mean RD"
	TARGET_FEATURE_NAMES = c(TARGET_FEATURE_NAMES, TARG_DEPTH)
	
	SAMP_DEPTH = "Sample mean RD"
	SAMPLE_FEATURE_NAMES = c(SAMPLE_FEATURE_NAMES, SAMP_DEPTH)
	
	if (PLOT_PC_CORRS) {
		for (job in JOBS) {
			PC = xhmm_data[[job]][["PC"]]
			PC_targs = colnames(PC)
			rm(PC)
		
			PC_LOADINGS = xhmm_data[[job]][["PC_LOADINGS"]]
			PC_LOADINGS_samps = colnames(PC_LOADINGS)
			rm(PC_LOADINGS)
		
			filteredRD = xhmm_data[[job]][["RD"]][PC_LOADINGS_samps, PC_targs]
		
			xhmm_data[[job]][[SAMP_DEPTH]] = rowMeans(filteredRD)
			xhmm_data[[job]][[TARG_DEPTH]] = colMeans(filteredRD)
		}
		
		
		
		PC.corr = list()
		for (job in JOBS) {
			PC = xhmm_data[[job]][["PC"]]
			PC_targs = colnames(PC)
			PC_t = t(PC)
			rm(PC)
		
			vals = matrix(0, nrow=length(PC_targs), ncol=length(TARGET_FEATURE_NAMES), dimnames=list(PC_targs, TARGET_FEATURE_NAMES))
			for (targVal in TARGET_FEATURE_NAMES) {
				vals[, targVal] = xhmm_data[[job]][[targVal]][PC_targs]
			}
			PC.corr[[job]] = cor(PC_t, vals, use="na.or.complete")
		}
		
		PC_LOADINGS.corr = list()
		for (job in JOBS) {
			PC_LOADINGS = xhmm_data[[job]][["PC_LOADINGS"]]
			# Get rid of the last principal component which has no standard deviation:
			PC_LOADINGS = PC_LOADINGS[1:(nrow(PC_LOADINGS)-1), ]
		
			PC_LOADINGS_samps = colnames(PC_LOADINGS)
			PC_LOADINGS_t = t(PC_LOADINGS)
			rm(PC_LOADINGS)
		
			vals = matrix(0, nrow=length(PC_LOADINGS_samps), ncol=length(SAMPLE_FEATURE_NAMES), dimnames=list(PC_LOADINGS_samps, SAMPLE_FEATURE_NAMES))
			for (sampVal in SAMPLE_FEATURE_NAMES) {
				vals[, sampVal] = xhmm_data[[job]][[sampVal]][PC_LOADINGS_samps]
			}
		
			# Add back in 0 correlation for last PC's loadings:
			PC_LOADINGS.corr[[job]] = rbind(cor(PC_LOADINGS_t, vals, use="na.or.complete"), rep(0, length(SAMPLE_FEATURE_NAMES)))
		}
		
		absCorrWithPC = list()
		
		pdf(paste(PLOT_PATH, "/PC_correlations.pdf", sep=""))
		PT_EXPANSION_FACTOR = 1.3
		for (job in JOBS) {
			numPC = nrow(PC.corr[[job]])
			NUM_PC_REMOVED = xhmm_data[[job]][["NUM_PC_REMOVED"]]
		
			absCorrWithPC[[job]] = abs(cbind(PC.corr[[job]], PC_LOADINGS.corr[[job]]))
		
			numTargFeatures = ncol(PC.corr[[job]])
			targFeatureColors = rainbow(numTargFeatures)
		
			numSampFeatures = ncol(PC_LOADINGS.corr[[job]])
			sampFeatureColors = rainbow(numSampFeatures)
		
			# 20 = small solid circles
			# 18 = small solid diamonds
			pch = c(rep(20, numTargFeatures), rep(18, numSampFeatures))
			colors = c(targFeatureColors, sampFeatureColors)
			featureNames = c(TARGET_FEATURE_NAMES, SAMPLE_FEATURE_NAMES)
		
			ylim = c(0,1)
		
			numPCsToPlot = c(numPC, 25, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500)
			for (plotPC in numPCsToPlot) {
				if (plotPC > numPC) {
					next
				}
				indsToPlot = 1:plotPC
				corrsToPlot = absCorrWithPC[[job]][indsToPlot, ]
				matplot(indsToPlot, corrsToPlot, ylim=ylim, pch=pch, type="p", col=colors, cex=PT_EXPANSION_FACTOR, xlab="Principal component", ylab="Absolute value of correlation to PC scores/loadings", main=paste(job, "", sep=""))
		
				plotY = seq(from=ylim[1], to=ylim[2], length.out=20)
				points(rep(NUM_PC_REMOVED, length(plotY)), plotY, col="black", type="p", pch="|")
		
				legend("top", featureNames, pch=pch, col=colors, pt.cex=PT_EXPANSION_FACTOR)
			}
		}
		dev.off()
		
		
		pdf(paste(PLOT_PATH, "/PC_factors.pdf", sep=""))

		oma = par("oma")
		# Add a space of 10 to the bottom outer margin to leave room for text on x-axis:
		oma[1] = oma[1] + 10
		par(oma=oma)

		for (job in JOBS) {
			NUM_PC_REMOVED = xhmm_data[[job]][["NUM_PC_REMOVED"]]
			absCorrWithRemovedPC = absCorrWithPC[[job]][1:NUM_PC_REMOVED, ]
		
			title = paste(job, ": correlation with PC", sep="")
			heatmap.2(absCorrWithRemovedPC, Rowv = NA, Colv = NA, dendrogram = "none", scale="none", revC = FALSE, density.info="none", trace="none", xlab = "", main=title, symbreaks=FALSE, colsep=1:ncol(absCorrWithRemovedPC), sepcolor="white")
		}
		dev.off()
		
		
		pdf(paste(PLOT_PATH, "/PC_stddev.pdf", sep=""))
		for (job in JOBS) {
			title = job
		
			SDdata = xhmm_data[[job]][["PC_SD"]]
		        # Use useInds to prevent "underflow" on log scale [remove outlier of tiny variance]:
		        useInds = which(SDdata > 1)
		
			ylim = c(min(SDdata[useInds]), max(SDdata[useInds]))
		
			plot(useInds, SDdata[useInds], ylim=ylim, log="y", xlab="Principal component", ylab="Standard deviation in component", main=title)
		
			NUM_PC_REMOVED = xhmm_data[[job]][["NUM_PC_REMOVED"]]
			plotY = exp(seq(from=log(ylim[1]), to=log(ylim[2]), length.out=20))
			points(rep(NUM_PC_REMOVED, length(plotY)), plotY, col="black", type="p", pch="|")
		}
		dev.off()
		
		
		
		CORR_THRESH = 0.1
		CORR_NUM_DECIMAL_PLACES = 2
	
		system(paste("mkdir -p ", PLOT_PATH, "/PC", sep=""))
		for (job in JOBS) {
			PC = xhmm_data[[job]][["PC"]]
			numPC = nrow(PC)
			chr = targetsToChrBp1Bp2(colnames(PC))[["chr"]]
		
			plotNames = c()
			for (i in 1:numPC) {
				plotName = paste(PLOT_PATH, "/PC/PC.", job, "_", i, sep="")
				plotNames = c(plotNames, plotName)
				
				title = paste("PC ", i, ": ", sep="")
				largeCorrelationNames = colnames(absCorrWithPC[[job]])[ which(absCorrWithPC[[job]][i, ] >= CORR_THRESH) ]
				if (length(largeCorrelationNames) == 0) {
					title = paste(title, "no correlations >= ", CORR_THRESH, sep="")
				}
				else {
					largeCorrelations = round(absCorrWithPC[[job]][i, largeCorrelationNames], CORR_NUM_DECIMAL_PLACES)
					title = paste(title, paste(largeCorrelationNames, ": ", largeCorrelations, collapse=", ", sep=""), sep="")
				}
				
				plotAllChromosomeValues(chr, PC[i, ], "Projection to prinicipal component", title, plotName)
			}
		
			if (system("which convert", ignore.stdout=TRUE, ignore.stderr=TRUE) == 0) {
				system(paste("convert ", paste(paste(plotNames, ".png", sep=""), collapse=" "), " ", PLOT_PATH, "/PC.", job, ".pdf", sep=""))
			}
		}
	}
	
	
	
	
	if (PLOT_ALL_CNVS) {
		PLOT_DIR = paste(PLOT_PATH, "/plot_CNV", sep="")
		system(paste("mkdir -p ", PLOT_DIR, sep=""))
		for (job in JOBS) {
			XCNV_CALLS = xhmm_data[[job]][["XCNV_CALLS"]]
			if (!is.null(USE_XCNV_TO_PLOT) && job %in% names(USE_XCNV_TO_PLOT)) {
				xcnvFile = USE_XCNV_TO_PLOT[job]
				if (xcnvFile != "") {
					XCNV_CALLS = loadXCNVcalls(xcnvFile)
				}
				else {
					XCNV_CALLS = matrix(nrow=0, ncol=0)
				}
			}

			if (nrow(XCNV_CALLS) > 0) {
			for (cnvInd in 1:nrow(XCNV_CALLS)) {
				cnv = XCNV_CALLS[cnvInd, ]
				sample = as.character(cnv["SAMPLE"])
				interval = as.character(cnv["INTERVAL"])
				SQ = as.numeric(cnv["Q_SOME"])
				type = as.character(cnv["CNV"])
				numTargets = as.numeric(cnv["NUM_TARG"])

				MARK_SAMPLES = sample
				if (!is.null(INCLUDE_PEDIGREE_SAMPLES)) {
					includeRows = which(INCLUDE_PEDIGREE_SAMPLES[["pedigree"]][, "#ID"] == sample | INCLUDE_PEDIGREE_SAMPLES[["pedigree"]][, "PAT"] == sample | INCLUDE_PEDIGREE_SAMPLES[["pedigree"]][, "MAT"] == sample)
					MARK_SAMPLES = union(sample, as.character(unlist(INCLUDE_PEDIGREE_SAMPLES[["pedigree"]][includeRows, c("#ID", "PAT", "MAT")])))
				}
	
				plotName = paste(PLOT_DIR, "/sample_", sample, ".", interval, ".SQ_", SQ, sep="")
	
				chrom_start_stop = targetsToChrBp1Bp2(interval)
				MARK_INTERVALS=list()
				MARK_INTERVALS[[paste("XHMM ", type, " (", numTargets, " targets)", sep="")]] = list(as.numeric(chrom_start_stop$bp1), as.numeric(chrom_start_stop$bp2), col="black", cex=0.8)
	
				for (addTargs in c(0, NUM_ADD_TARGS)) {
					start_stop_inds = regionToStartStopInds(xhmm_data[[job]], interval, NUM_ADD_TARGS=addTargs)
	
					usePlotName = plotName
					if (addTargs == 0) {
						usePlotName = paste(usePlotName, ".exact", sep="")
					}
					plot_XHMM_targets(usePlotName, xhmm_data[[job]], xhmm_data[[job]][["TARGETS_TO_GENES"]], SAMPLE_FEATURES, SQ_THRESH, start_stop_inds[1], start_stop_inds[2], PLOT_ONLY_PNG=PLOT_ONLY_PNG, MARK_SAMPLES=MARK_SAMPLES, MARK_INTERVALS=MARK_INTERVALS, APPEND_REGION_NAME=FALSE)
				}
			}
			}
		}
	}

	if (PLOT_READ_DEPTHS) {
		pdf(paste(PLOT_PATH, "/per_target_sd.pdf", sep=""))
		for (job in JOBS) {
			pcaNormalizedRD = xhmm_data[[job]][["PCA_NORMALIZED"]]
			sd.target = apply(pcaNormalizedRD, 2, sd)
			hist(sd.target, breaks=400, xlab="Per-target standard deviation of PCA-normalized read depth (coverage)", ylab="Frequency", main=paste(job, sep=""))
		}
		dev.off()
	}

	return(xhmm_data)
}
