MAIT_identify_metabolites = function(dataset, metadata.variable, xSet = NULL, data.folder = NULL, features = NULL, 
									 mass.tolerance = 0.5){

	metadata.var = dataset$metadata[,metadata.variable]
	mait.object = sampleProcessing_modified(dataDir = data.folder, metadata = metadata.var, xSet = xSet, project = "MAIT")
	mait.annotation = MAIT::peakAnnotation(MAIT.object = mait.object, corrWithSamp = 0.7, corrBetSamp = 0.75, 
                            perfwhm = 0.6)
	mait.sig = MAIT::spectralSigFeatures(MAIT.object = mait.annotation, pvalue = 0.05, p.adj = "none",
                               scale = FALSE, printCSVfile = F)
    
	scoresTable = MAIT::getScoresTable(MAIT.object = mait.sig, getExtendedTable = T)

	if (!is.null(features) && length(features) == 1 && features == "all"){
		mait.sig@FeatureData@featureSigID = 1:length(scoresTable$extendedTable$mz)
	} else if (!is.null(features) && length(features) > 1){
		mzs = scoresTable$extendedTable$mz
		indexes = which(mzs > (features - mass.tolerance) & mzs < (features + mass.tolerance))
		mait.sig@FeatureData@featureSigID = indexes
    }
	
	mait.bio = MAIT::Biotransformations(MAIT.object = mait.sig, adductAnnotation = T, peakPrecision = 0.005)
	mait.identify = MAIT::identifyMetabolites(MAIT.object = mait.sig, peakTolerance = 0.005)
	mait.identify
}

sampleProcessing_modified = function (dataDir = NULL, metadata = NULL, xSet = NULL, snThres = 2, Sigma = 5/2.3548, mzSlices = 0.3, 
    retcorrMethod = "loess", groupMethod = "density", bwGroup = 3, 
    mzWidGroup = 0.25, filterMethod = "matchedFilter", rtStep = 0.03, 
    nSlaves = 0, minfrac = 0.5, minsamp = 1, peakwidth = c(5, 
        20), project = NULL, ppm = 10, family = c("gaussian", 
        "symmetric"), span = 0.2, fwhm = 30) 
{
    if (is.null(dataDir)) {
        stop("No input directory was given")
    }
    if (is.null(project)) {
        stop("No project name was included")
    }
    parameters <- list(dataDir, snThres, Sigma, mzSlices, retcorrMethod, 
        groupMethod, bwGroup, mzWidGroup, filterMethod, rtStep, 
        nSlaves, project, ppm, minfrac, fwhm, family, span, peakwidth)
    names(parameters) <- c("dataDir", "snThres", "Sigma", "mzSlices", 
        "retcorrMethod", "groupMethod", "bwGroup", "mzWidGroup", 
        "filterMethod", "rtStep", "nSlaves", "project", "ppm", 
        "minfrac", "fwhm", "family", "span", "centWave peakwidth")
    MAIT.object <- new("MAIT")
    MAIT.object@RawData@parameters@sampleProcessing <- parameters		
    fileList <- list.files(dataDir, recursive = T, full.names = TRUE)
    #writeParameterTable(parameters(MAIT.object), folder = resultsPath(MAIT.object))

    if (is.null(metadata)){
		class <- c("Class1", "Class2")
		classNum <- vector(length = length(class))
		firstClass = length(fileList)%/%2
		secondClass = length(fileList) - firstClass
		classNum[1] <- firstClass
		classNum[2] <- secondClass 
		classes <- rep(class, classNum)
	} else {
		class = levels(metadata)
		classNum <- vector(length = length(table(class)))
		for (i in 1:length(classNum)){
			classNum[i] = as.numeric(table(metadata)[i])
		}
		classes <- as.character(metadata)
	}
	
    if (length(list.files(dataDir)) == 1) {
        warning("Warning: Input data only has one class!")
    }
    if (is.null(project)) {
        warning("Warning: Project name is empty!")
    }
    if (!is.null(project)) {
        resultsPath <- paste("Results", project, sep = "_")
        dir.create(resultsPath)
    }
    else {
        resultsPath <- "Results"
        dir.create(resultsPath)
    }
    
    if (!is.null(xSet)){
		#xSet@phenoData$class = classes
		snames = rownames(xSet@phenoData)
		xcms::phenoData(xSet) = classes
		rownames(xcms::phenoData(xSet)) <- snames
		fPeaks = list(xSet)
    } else {
		if (filterMethod == "matchedFilter") {
			peaks <- xcms::xcmsSet(files = fileList, snthresh = snThres, 
				method = filterMethod, sigma = Sigma, max = 3, step = rtStep, 
				mzdiff = mzSlices, sclass = classes, nSlaves = nSlaves, 
				fwhm = fwhm)
		} 	
		if (filterMethod == "centWave") {
			peaks <- xcms::xcmsSet(files = fileList, snthresh = snThres, 
				method = filterMethod, ppm = ppm, mzdiff = mzSlices, 
				sclass = classes, nSlaves = nSlaves, peakwidth = peakwidth)
		}
		
		cat("Peak detection done", fill = TRUE)
		groups <- xcms::group(peaks, method = groupMethod, bw = bwGroup, 
			mzwid = mzWidGroup, max = 50, minfrac = minfrac, minsamp = minsamp)
		if (retcorrMethod != "none") {
			retcorr_groups <- xcms::retcor(groups, method = retcorrMethod, family = family, span = span)
			cat("Retention time correction done", fill = TRUE)
			groups <- xcms::group(retcorr_groups, method = groupMethod, 
				bw = bwGroup, mzwid = mzWidGroup, max = 50)
			cat("Peak grouping after samples done", fill = TRUE)
		}
		else {
			cat("Skipping retention time correction...", fill = TRUE)
		}
		fPeaks <- xcms::fillPeaks(groups)
		cat("Missing Peak integration done", fill = TRUE)
		fPeaks <- list(fPeaks)
	}


    names(fPeaks) <- "xcmsSet"
    MAIT.object@RawData@data <- fPeaks
    MAIT.object@PhenoData@classes <- class
    MAIT.object@PhenoData@classNum <- classNum
    MAIT.object@PhenoData@resultsPath <- resultsPath
    return(MAIT.object)
}
