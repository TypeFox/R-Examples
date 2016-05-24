# 
# mzmatch_identify_metabolites = function(ionisation="detect", data.folder = NULL, xSet = NULL, backend = "Ramp",
# 										adducts = "M+H,M+ACN+Na,M+Na,M+K,M+ACN+H", databases = NULL){
# 										
# 	mzmatch.R::mzmatch.init(version.1=FALSE)
# 	
# 	if (is.null(xSet)){
# 		files = list.files(data.folder, recursive=T, full.names=TRUE)
# 		xSet <- xcms::xcmsSet(files, method='centWave', ppm=2, peakwidth=c(5,100), 
#                   snthresh=5, prefilter=c(3,1000), integrate=1, mzdiff=0.001, 
#                   verbose.columns=TRUE, fitgauss=FALSE, nSlaves=4)
# 	}
# 	
# 	peakml.files = gsub("(.+)/([^/]+\\..+$)","\\2", xSet@filepaths)
# 	peakml.files = gsub("\\..+$", "\\.peakml", peakml.files)
# 	peakml.files = paste("./peakml/", peakml.files, sep="")
# 	dir.create("peakml")
# 	mzmatch.R::PeakML.xcms.write.SingleMeasurement (xset=xSet,filename=peakml.files,
#                                        ionisation=ionisation,addscans=2,
#                                        writeRejected=FALSE,ApodisationFilter=TRUE)
#     
# 	mzmatch.ipeak.Combine_modified(i=paste(peakml.files,collapse=","),v=T,rtwindow=30,
# 						  o="./mzmatch/combined.peakml",combination="set",ppm=5)
# 	mzmatch.ipeak.filter.NoiseFilter_modified (i="./mzmatch/combined.peakml",o="./mzmatch/combined_noisef.peakml",
# 									  v=T,codadw=0.8)
# 	mzmatch.ipeak.filter.SimpleFilter_modified(i="./mzmatch/combined_noisef.peakml", 
# 									  o="./mzmatch/combined_sfdet.peakml", mindetections=3)
# 	mzmatch.ipeak.filter.SimpleFilter_modified(i="./mzmatch/combined_sfdet.peakml", 
# 									  o="./mzmatch/combined_highintensity.peakml", 
# 									  minintensity=100000)
# 	PeakML.GapFiller_modified(filename = "./mzmatch/combined_highintensity.peakml", ionisation = ionisation, 
# 					 outputfile = "./mzmatch/highintensity_gapfilled.peakml", 
# 					 ppm=5, rtwin = 0, backend = backend)
# 	mzmatch.ipeak.sort.RelatedPeaks_modified (i="./mzmatch/highintensity_gapfilled.peakml",v=T,
# 									 o="./mzmatch/mzMatch_output.peakml",
# 									 basepeaks="./mzmatch/mzMatch_basepeaks.peakml",ppm=3,
# 									 rtwindow=6)
# 									 
# 	annot <- paste("relation.id,relation.ship,codadw,charge")
# 	mzmatch.ipeak.convert.ConvertToText_modified (i="./mzmatch/mzMatch_output.peakml",
# 										 o="./mzmatch/mzMATCHoutput.txt",v=T,annotations=annot)
# 
# 	if (is.null(databases)){
# 		DBS <- dir(paste(find.package("mzmatch.R"), "/dbs", sep=""),
# 			   full.names=TRUE)
# 		DBS <- paste(DBS,collapse='","')
# 		DBS <- paste('"',DBS,'"',sep="")
# 	} else {
# 		DBS = databases
# 	}
# 
# 
# 	mzmatch.ipeak.util.Identify_modified(i="./mzmatch/mzMatch_output.peakml", v=T,
# 								o="./mzmatch/metabolites.peakml", ppm=3, 
# 								databases=DBS, adducts = adducts)
# 
# 
# 	mzmatch.ipeak.convert.ConvertToText_modified (
# 	  i="./mzmatch/metabolites.peakml",
# 	  o= "./mzmatch/metabolites.txt", databases=DBS,
# 	  annotations="identification,moleculeName,ppm,adduct,relation.ship")
#     
# 	metabolites = mzmatch.R::PeakML.Read("./mzmatch/metabolites.peakml")
# 	metabolites
# }
# 
# 
# PeakML.GapFiller_modified = function (filename, ionisation = "detect", Rawpath = NULL, outputfile, 
#     ppm = 0, rtwin = 0, nSlaves = 1, fillAll = FALSE, backend = "Ramp") 
# {
#     pos = 1
#     version.1 <- get("version.1", envir = as.environment(pos))
#     FillinPeaks <- function(peaknum) {
#         whichpeakset <- numchromsexpected[fillinnums[peaknum], 
#             1]
#         subtable <- PeakMLdata$peakDataMtx[PeakMLdata$peakDataMtx[, 
#             10] == whichpeakset, ]
#         subtable <- rbind(subtable, NULL)
#         rt_start <- min(subtable[, 5]) - rtwin
#         rt_finis <- max(subtable[, 6]) + rtwin
#         if (rt_finis > max(correctedRT)) {
#             rt_finis <- max(correctedRT)
#         }
#         if (rt_start > max(correctedRT)) {
#             rt_start <- max(correctedRT)
#         }
#         if (rt_finis < min(correctedRT)) {
#             rt_finis <- min(correctedRT)
#         }
#         if (rt_start < min(correctedRT)) {
#             rt_start <- min(correctedRT)
#         }
#         mz_start <- min(subtable[, 2])
#         mz_finis <- max(subtable[, 3])
#         mz_start <- mz_start - (mz_start * ppm/10^6)
#         mz_finis <- mz_finis + (mz_finis * ppm/10^6)
#         scan_start <- which(correctedRT >= rt_start)[1] - 1
#         if (scan_start == 0) 
#             scan_start = 1
#         scan_finis <- which(correctedRT >= rt_finis)[1]
#         C <- try(mzmatch.R::PeakML.Methods.getRawMat(allRawPeaks, scan_start, 
#             scan_finis, mz_start, mz_finis, correctedRT, uncorrectedRT))
#         if (class(C) == "try-error") {
#             C <- c(1, 1, 1, 1, 1)
#         }
#         C <- rbind(C, NULL)
#         if (nrow(C) <= 3 | length(unique(C[, 5])) <= 3) {
#             scanids <- c(-1, -1, -1)
#             retentiontimes <- c(-1, -1, -1)
#             masses <- c(-1, -1, -1)
#             intensities <- c(-1, -1, -1)
#         }
#         else {
#             scanids <- C[, 3]
#             retentiontimes <- C[, 2]
#             masses <- C[, 4]
#             intensities <- C[, 5]
#         }
#         OUT <- rbind(masses, intensities, retentiontimes, scanids - 
#             1)
#     }
#     st <- system.time(PeakMLdata <- mzmatch.R::PeakML.Read(filename, ionisation, 
#         Rawpath))
#     ionisation <- PeakMLdata$massCorrection[[2]]
#     massCorrection <- PeakMLdata$massCorrection[[1]]
#     samplenames <- PeakMLdata$sampleNames
#     rawdatafullpaths <- PeakMLdata$rawDataFullPaths
#     if (is.null(rawdatafullpaths)) {
#         cat("Some of the raw data files are not accessible, we will not be able to fill in missing peaks. Please set \"Rawpath\" argument with location where files can be located\n")
#         stop()
#     }
#     numchromsexpected <- unlist(lapply(1:max(PeakMLdata$peakDataMtx[, 
#         10]), function(x) rep(x, length(samplenames))))
#     numchromsexpected <- cbind(numchromsexpected, NA, NA, NA)
#     for (setnum in 1:max(PeakMLdata$peakDataMtx[, 10])) {
#         inset <- c(1:length(samplenames))
#         rownums <- which(PeakMLdata$peakDataMtx[, 10] == setnum)
#         hit <- PeakMLdata$peakDataMtx[rownums, 9]
#         oneEqualsSetnum = which(numchromsexpected[, 1] == setnum)
#         numchromsexpected[oneEqualsSetnum, 2] <- as.numeric(inset %in% 
#             hit)
#         missed <- which(inset %in% hit == FALSE)
#         if (length(missed) > 0) {
#             detectedpeaks <- c(rep(1, length(hit)), rep(0, length(missed)))
#             hit <- append(hit, missed)
#             rownums <- append(rownums, rep(0, length(missed)))
#         }
#         else {
#             detectedpeaks <- rep(1, length(hit))
#         }
#         numchromsexpected[oneEqualsSetnum, 3] <- hit
#         numchromsexpected[oneEqualsSetnum, 2] <- detectedpeaks
#         numchromsexpected[oneEqualsSetnum, 4] <- rownums
#     }
#     colnames(numchromsexpected) <- NULL
#     chromslist <- vector("list", nrow(numchromsexpected))
#     if (fillAll == TRUE) {
#         detectedchromatograms <- numchromsexpected[, 2]
#         numchromsexpected[, 2] <- 0
#     }
#     notdetected <- which(numchromsexpected[, 2] == 0)
#     whichfiles <- numchromsexpected[notdetected, 3]
#     samplenums <- unique(whichfiles)
#     if (length(samplenums != 0)) {
#         for (filenum in 1:length(samplenums)) {
#             samplefile <- samplenums[filenum]
#             cat("Working on file: ", rawdatafullpaths[samplefile], 
#                 "\n")
#             rawfile <- mzR::openMSfile(rawdatafullpaths[samplefile], 
#                 verbose = FALSE, backend = backend)
#             allRawPeaks <- mzR::peaks(rawfile)
#             correctedRT <- as.numeric(PeakMLdata$correctedRTList[[samplefile]])
#             uncorrectedRT <- mzR::header(rawfile)$retentionTime
#             if (all(correctedRT == uncorrectedRT)) {
#                 rtCorrection <- FALSE
#             }
#             else {
#                 rtCorrection <- TRUE
#             }
#             fillinnums <- notdetected[whichfiles == samplefile]
#             isSnow <- FALSE
#             if (nSlaves > 1) {
#                 HIT <- grep("parallel", installed.packages()[, 1])
#                 if (!is.null(HIT)) {
#                   isSnow == TRUE
#                 }
#                 else {
#                   cat("Pleae install package parallel to use multiple processors. \n We will continue with a single processor for the time being.", 
#                     "\\n")
#                 }
#             }
#             if (isSnow == TRUE) {
#                 if (filenum == 1) {
#                   cat("Package parallel loaded.", "\n")
#                 }
# 		rawMat = NULL
#                 cl <- parallel::makeCluster(nSlaves)
#                 assign("rtwin", rtwin, envir = as.environment(pos))
#                 assign("rawfile", rawfile, envir = as.environment(pos))
#                 assign("numchromsexpected", numchromsexpected, 
#                   envir = as.environment(pos))
#                 assign("fillinnums", fillinnums, envir = as.environment(pos))
#                 assign("PeakMLdata$peakDataMtx", PeakMLdata$peakDataMtx, 
#                   envir = as.environment(pos))
#                 assign("ppm", ppm, envir = as.environment(pos))
#                 assign("FillinPeaks", FillinPeaks, ppm, envir = as.environment(pos))
#                 assign("rawMat", rawMat, envir = as.environment(pos))
#                 parallel::clusterExport(cl, list = c("rtwin", "rawfile", 
#                   "numchromsexpected", "fillinnums", "PeakMLdata$peakDataMtx", 
#                   "ppm", "FillinPeaks", "rawMat"))
#                 filledlist <- parallel::parLapply(cl, c(1:length(fillinnums)), 
#                   FillinPeaks)
#                 parallel::stopCluster(cl)
#             }
#             else {
#                 filledlist <- lapply(1:length(fillinnums), FillinPeaks)
#             }
#             for (i in 1:length(filledlist)) {
#                 chromslist[[fillinnums[i]]] <- filledlist[[i]]
#             }
#             mzR::close(rawfile)
#             rm(filledlist)
#         }
#     }
#     project <- rJava::.jnew("peakml/util/rjava/Project", samplenames, 
#         rawdatafullpaths, as.character(PeakMLdata$phenoData))
#     rJava::.jcall(project, returnSig = "V", method = "addHeaderAnnotation", 
#         as.character("peakproc"), as.character("XCMS_Gapfilled"))
#     for (measurementid in 1:length(samplenames)) {
#         for (scannum in 1:length(PeakMLdata$correctedRTList[[measurementid]])) {
#             rJava::.jcall(project, returnSig = "V", method = "addScanInfo", 
#                 as.integer(measurementid - 1), as.numeric(PeakMLdata$correctedRTList[[measurementid]][scannum]), 
#                 as.character(ionisation))
#             rJava::.jcall(project, returnSig = "V", method = "addScanAnnotation", 
#                 as.integer(measurementid - 1), as.integer(scannum - 
#                   1), as.character("RT_raw"), as.character(PeakMLdata$rawRTList[[measurementid]][scannum]))
#         }
#     }
#     for (i in 1:length(chromslist)) {
#         if (numchromsexpected[i, 2] == 0) {
#             chrom <- chromslist[[i]]
#         }
#         else {
#             ind <- numchromsexpected[i, 4]
#             chrom <- PeakMLdata$chromDataList[[ind]]
#         }
#         rJava::.jcall(project, returnSig = "V", method = "addMassChromatogram", 
#             as.integer(numchromsexpected[i, 3] - 1), as.integer(chrom[4, 
#                 ]), as.numeric(chrom[3, ]), as.numeric(chrom[1, 
#                 ]), as.numeric(chrom[2, ]), as.character(ionisation))
#     }
#     setindexes <- vector("list", length(unique(numchromsexpected[, 
#         1])))
#     for (indexnumber in 1:length(setindexes)) {
#         setindexes[[indexnumber]] <- which(numchromsexpected[, 
#             1] == indexnumber)
#     }
#     for (ind in 1:length(setindexes)) {
#         rJava::.jcall(project, returnSig = "V", method = "addPeakSet", 
#             as.integer(setindexes[[ind]] - 1))
#     }
#     if (!is.null(PeakMLdata$GroupAnnotations)) {
#         mzmatch.R::PeakML.Methods.writeGroupAnnotations(project, PeakMLdata$GroupAnnotations)
#     }
#     rJava::.jcall(project, returnSig = "V", method = "write", outputfile)
# }
# 
# 
# mzmatch.ipeak.util.Identify_modified = function (JHeapSize = 1425, i = NULL, o = NULL, ppm = NULL, databases = NULL, 
#     minrt = NULL, maxrt = NULL, rtwindow = NULL, rtwindowrelative = NULL, 
#     massOverride = NULL, polarity = NULL, adducts = NULL, h = NULL, 
#     v = NULL) 
# {
#     pos = 1
#     version.1 <- get("version.1", envir = as.environment(pos))
#     java <- "java -da -dsa -Xmn1g -Xss228k -XX:+UseParallelGC -XX:ParallelGCThreads=10"
#     JHeapSize <- paste(JHeapSize, "m", sep = "")
#     java <- paste(java, " -Xms", JHeapSize, " -Xmx", JHeapSize, 
#         " -cp", sep = "")
#     if (version.1 == TRUE) {
#         mzmatch <- paste(java, ' "', find.package("mzmatch.R"), 
#             "/java/mzmatch.jar", '"', sep = "")
#     }
#     else {
#         mzmatch <- paste(java, ' "', find.package("mzmatch.R"), 
#             "/java/mzmatch_2.0.jar", '"', sep = "")
#     }
#     tool <- paste(mzmatch, "mzmatch.ipeak.util.Identify")
#     if (!is.null(i)) 
#         tool <- paste(tool, "-i", i)
#     if (!is.null(o)) 
#         tool <- paste(tool, "-o", o)
#     if (!is.null(ppm)) 
#         tool <- paste(tool, "-ppm", ppm)
#     if (!is.null(minrt)) 
#         tool <- paste(tool, "-minrt", minrt)
#     if (!is.null(maxrt)) 
#         tool <- paste(tool, "-maxrt", maxrt)
#     if (!is.null(rtwindow)) 
#         tool <- paste(tool, "-rtwindow", rtwindow)
#     if (!is.null(rtwindowrelative) && rtwindowrelative == T && 
#         !is.null(rtwindow)) 
#         tool <- paste(tool, "-rtwindowrelative true")
#     if (!is.null(massOverride)) 
#         tool <- paste(tool, "-massOverride", massOverride)
#     if (!is.null(polarity)) 
#         tool <- paste(tool, "-polarity", polarity)
#     if (!is.null(adducts)) 
#         tool <- paste(tool, "-adducts", adducts)
# 	if (!is.null(databases)) 
#         tool <- paste(tool, "-databases", databases)
#     if (!is.null(h) && h == T) 
#         tool <- paste(tool, "-h")
#     if (!is.null(v) && v == T) 
#         tool <- paste(tool, "-v")
#     cat(tool)
# 	print(tool)
#     system(tool)
# }
# 
# mzmatch.ipeak.convert.ConvertToText_modified = function (JHeapSize = 1425, i = NULL, o = NULL, databases = NULL, 
#     annotations = NULL, h = NULL, v = NULL) 
# {
#     pos = 1
#     version.1 <- get("version.1", envir = as.environment(pos))
#     java <- "java -da -dsa -Xmn1g -Xss228k -XX:+UseParallelGC -XX:ParallelGCThreads=10"
#     JHeapSize <- paste(JHeapSize, "m", sep = "")
#     java <- paste(java, " -Xms", JHeapSize, " -Xmx", JHeapSize, 
#         " -cp", sep = "")
#     if (version.1 == TRUE) {
#         mzmatch <- paste(java, ' "', find.package("mzmatch.R"), 
#             "/java/mzmatch.jar", '"', sep = "")
#     }
#     else {
#         mzmatch <- paste(java, ' "', find.package("mzmatch.R"), 
#             "/java/mzmatch_2.0.jar", '"', sep = "")
#     }
#     tool <- paste(mzmatch, "mzmatch.ipeak.convert.ConvertToText")
#     if (!is.null(i)) 
#         tool <- paste(tool, "-i", i)
#     if (!is.null(o)) 
#         tool <- paste(tool, "-o", o)
#     if (!is.null(databases)) 
#         tool <- paste(tool, "-databases", databases)
#     if (!is.null(annotations)) 
#         tool <- paste(tool, "-annotations", annotations)
#     if (!is.null(h) && h == T) 
#         tool <- paste(tool, "-h")
#     if (!is.null(v) && v == T) 
#         tool <- paste(tool, "-v")
#     system(tool)
# }
# 
# mzmatch.ipeak.sort.RelatedPeaks_modified = function (JHeapSize = 1425, i = NULL, o = NULL, basepeaks = NULL, 
#     ppm = NULL, rtwindow = NULL, minrt = NULL, h = NULL, v = NULL) 
# {
#     pos = 1
#     version.1 <- get("version.1", envir = as.environment(pos))
#     java <- "java -da -dsa -Xmn1g -Xss228k -XX:+UseParallelGC -XX:ParallelGCThreads=10"
#     JHeapSize <- paste(JHeapSize, "m", sep = "")
#     java <- paste(java, " -Xms", JHeapSize, " -Xmx", JHeapSize, 
#         " -cp", sep = "")
#     if (version.1 == TRUE) {
#         mzmatch <- paste(java, ' "', find.package("mzmatch.R"), 
#             "/java/mzmatch.jar", '"', sep = "")
#     }
#     else {
#         mzmatch <- paste(java, ' "', find.package("mzmatch.R"), 
#             "/java/mzmatch_2.0.jar", '"', sep = "")
#     }
#     tool <- paste(mzmatch, "mzmatch.ipeak.sort.RelatedPeaks")
#     if (!is.null(i)) 
#         tool <- paste(tool, "-i", i)
#     if (!is.null(o)) 
#         tool <- paste(tool, "-o", o)
#     if (!is.null(basepeaks)) 
#         tool <- paste(tool, "-basepeaks", basepeaks)
#     if (!is.null(ppm)) 
#         tool <- paste(tool, "-ppm", ppm)
#     if (!is.null(rtwindow)) 
#         tool <- paste(tool, "-rtwindow", rtwindow)
#     if (!is.null(minrt)) 
#         tool <- paste(tool, "-minrt", minrt)
#     if (!is.null(h) && h == T) 
#         tool <- paste(tool, "-h")
#     if (!is.null(v) && v == T) 
#         tool <- paste(tool, "-v")
#     system(tool)
# }
# 
# mzmatch.ipeak.filter.SimpleFilter_modified = function (JHeapSize = 1425, i = NULL, o = NULL, rejected = NULL, 
#     databases = NULL, ppm = NULL, n = NULL, offset = NULL, mindetections = NULL, 
#     minscanid = NULL, maxscanid = NULL, minretentiontime = NULL, 
#     maxretentiontime = NULL, minmass = NULL, maxmass = NULL, 
#     minintensity = NULL, maxintensity = NULL, annotations = NULL, 
#     h = NULL, v = NULL) 
# {
#     pos = 1
#     version.1 <- get("version.1", envir = as.environment(pos))
#     java <- "java -da -dsa -Xmn1g -Xss228k -XX:+UseParallelGC -XX:ParallelGCThreads=10"
#     JHeapSize <- paste(JHeapSize, "m", sep = "")
#     java <- paste(java, " -Xms", JHeapSize, " -Xmx", JHeapSize, 
#         " -cp", sep = "")
#     if (version.1 == TRUE) {
#         mzmatch <- paste(java, ' "', find.package("mzmatch.R"), 
#             "/java/mzmatch.jar", '"', sep = "")
#     }
#     else {
#         mzmatch <- paste(java, ' "', find.package("mzmatch.R"), 
#             "/java/mzmatch_2.0.jar", '"', sep = "")
#     }
#     tool <- paste(mzmatch, "mzmatch.ipeak.filter.SimpleFilter")
#     if (!is.null(i)) 
#         tool <- paste(tool, "-i", i)
#     if (!is.null(o)) 
#         tool <- paste(tool, "-o", o)
#     if (!is.null(rejected)) 
#         tool <- paste(tool, "-rejected", rejected)
#     if (!is.null(databases)) 
#         tool <- paste(tool, "-databases", databases)
#     if (!is.null(ppm)) 
#         tool <- paste(tool, "-ppm", ppm)
#     if (!is.null(n)) 
#         tool <- paste(tool, "-n", n)
#     if (!is.null(offset)) 
#         tool <- paste(tool, "-offset", offset)
#     if (!is.null(mindetections)) 
#         tool <- paste(tool, "-mindetections", mindetections)
#     if (!is.null(minscanid)) 
#         tool <- paste(tool, "-minscanid", minscanid)
#     if (!is.null(maxscanid)) 
#         tool <- paste(tool, "-maxscanid", maxscanid)
#     if (!is.null(minretentiontime)) 
#         tool <- paste(tool, "-minretentiontime", minretentiontime)
#     if (!is.null(maxretentiontime)) 
#         tool <- paste(tool, "-maxretentiontime", maxretentiontime)
#     if (!is.null(minmass)) 
#         tool <- paste(tool, "-minmass", minmass)
#     if (!is.null(maxmass)) 
#         tool <- paste(tool, "-maxmass", maxmass)
#     if (!is.null(minintensity)) 
#         tool <- paste(tool, "-minintensity", minintensity)
#     if (!is.null(maxintensity)) 
#         tool <- paste(tool, "-maxintensity", maxintensity)
#     if (!is.null(annotations)) 
#         tool <- paste(tool, "-annotations", annotations)
#     if (!is.null(h) && h == T) 
#         tool <- paste(tool, "-h")
#     if (!is.null(v) && v == T) 
#         tool <- paste(tool, "-v")
#     system(tool)
# }
# 
# mzmatch.ipeak.filter.NoiseFilter_modified = function (JHeapSize = 1425, i = NULL, o = NULL, rejected = NULL, 
#     codadw = NULL, h = NULL, v = NULL) 
# {
#     pos = 1
#     version.1 <- get("version.1", envir = as.environment(pos))
#     java <- "java -da -dsa -Xmn1g -Xss228k -XX:+UseParallelGC -XX:ParallelGCThreads=10"
#     JHeapSize <- paste(JHeapSize, "m", sep = "")
#     java <- paste(java, " -Xms", JHeapSize, " -Xmx", JHeapSize, 
#         " -cp", sep = "")
#     if (version.1 == TRUE) {
#         mzmatch <- paste(java, ' "', find.package("mzmatch.R"), 
#             "/java/mzmatch.jar", '"',sep = "")
#     }
#     else {
#         mzmatch <- paste(java, ' "', find.package("mzmatch.R"), 
#             "/java/mzmatch_2.0.jar", '"', sep = "")
#     }
#     tool <- paste(mzmatch, "mzmatch.ipeak.filter.NoiseFilter")
#     if (!is.null(i)) 
#         tool <- paste(tool, "-i", i)
#     if (!is.null(o)) 
#         tool <- paste(tool, "-o", o)
#     if (!is.null(rejected)) 
#         tool <- paste(tool, "-rejected", rejected)
#     if (!is.null(codadw)) 
#         tool <- paste(tool, "-codadw", codadw)
#     if (!is.null(h) && h == T) 
#         tool <- paste(tool, "-h")
#     if (!is.null(v) && v == T) 
#         tool <- paste(tool, "-v")
#     system(tool)
# }
# 
# mzmatch.ipeak.Combine_modified= function (JHeapSize = 1425, i = NULL, o = NULL, label = NULL, 
#     labels = NULL, ppm = NULL, rtwindow = NULL, combination = NULL, 
#     h = NULL, v = NULL, sampleList = NULL, nSlaves = 1, outputfolder = "combined") 
# {
#     pos = 1
#     version.1 <- get("version.1", envir = as.environment(pos))
#     java <- "java -da -dsa -Xmn1g -Xss228k -XX:+UseParallelGC -XX:ParallelGCThreads=10"
#     JHeapSize <- paste(JHeapSize, "m", sep = "")
#     java <- paste(java, " -Xms", JHeapSize, " -Xmx", JHeapSize, 
#         " -cp", sep = "")
#     if (!is.null(i)) {
#         if (version.1 == TRUE) {
#             mzmatch <- paste(java, ' "', find.package("mzmatch.R"), 
#                 "/java/mzmatch.jar", '"', sep = "")
#         }
#         else {
#             mzmatch <- paste(java, ' "', find.package("mzmatch.R"), 
#                 "/java/mzmatch_2.0.jar", '"',sep = "")
#         }
#         tool <- paste(mzmatch, "mzmatch.ipeak.Combine")
#         if (!is.null(i)) 
#             tool <- paste(tool, "-i", i)
#         if (!is.null(o)) 
#             tool <- paste(tool, "-o", o)
#         if (!is.null(label)) 
#             tool <- paste(tool, "-label", label)
#         if (!is.null(labels)) 
#             tool <- paste(tool, "-labels", labels)
#         if (!is.null(ppm)) 
#             tool <- paste(tool, "-ppm", ppm)
#         if (!is.null(rtwindow)) 
#             tool <- paste(tool, "-rtwindow", rtwindow)
#         if (!is.null(combination)) 
#             tool <- paste(tool, "-combination", combination)
#         if (!is.null(h) && h == T) 
#             tool <- paste(tool, "-h")
#         if (!is.null(v) && v == T) 
#             tool <- paste(tool, "-v")
#         system(tool)
#     }
#     else {
#         MainClasses <- levels(as.factor(sampleList$sampleClass))
#         if (!file.exists(outputfolder)) {
#             dir.create(outputfolder)
#         }
#         combineFunction <- function(fnum) {
#             FILESf <- sampleList$outputfilenames[which(sampleList$sampleClass == 
#                 MainClasses[fnum])]
#             OUTPUTf <- paste(outputfolder, "/", MainClasses[fnum], 
#                 ".peakml", sep = "")
#             mzmatch <- paste(java, " ", find.package("mzmatch.R"), 
#                 "/java/mzmatch_2.0.jar", sep = "")
#             tool <- paste(mzmatch, "mzmatch.ipeak.Combine")
#             i = paste(FILESf, collapse = ",")
#             tool <- paste(tool, "-i", i)
#             o = OUTPUTf
#             tool <- paste(tool, "-o", o)
#             label = MainClasses[fnum]
#             tool <- paste(tool, "-label", label)
#             if (!is.null(labels)) 
#                 tool <- paste(tool, "-labels", labels)
#             if (!is.null(ppm)) 
#                 tool <- paste(tool, "-ppm", ppm)
#             if (!is.null(rtwindow)) 
#                 tool <- paste(tool, "-rtwindow", rtwindow)
#             if (!is.null(combination)) 
#                 tool <- paste(tool, "-combination", combination)
#             if (!is.null(h) && h == T) 
#                 tool <- paste(tool, "-h")
#             if (!is.null(v) && v == T) 
#                 tool <- paste(tool, "-v")
#             cat(label, "\n")
#             system(tool)
#         }
#         if (nSlaves > 1) {
#             cl <- parallel::makeCluster(nSlaves)
#             envname <- environment()
#             parallel::clusterExport(cl, list = c("combineFunction", "MainClasses", 
#                 "sampleList", "labels", "ppm", "rtwindow", "combination", 
#                 "java", "h", "v"), envir = envname)
#             system.time(parallel::clusterApply(cl, 1:length(MainClasses), 
#                 combineFunction))
#             parallel::stopCluster(cl)
#         }
#         else {
#             for (fnum in 1:length(MainClasses)) {
#                 combineFunction(fnum)
#             }
#         }
#     }
# }
