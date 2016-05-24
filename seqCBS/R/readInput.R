readInput <-
function(inputFilename, formatName="Chiang", sep="\t") {
	inputFiles = readListInputFile(inputFilename, sep)
	nNormalRun = length(inputFiles$normalFiles)
	nTumorRun = length(inputFiles$tumorFiles)
	normalList = vector("list", nNormalRun)
	tumorList = vector("list", nTumorRun)
	for(i in 1:nNormalRun) {
		normalList[[i]] = readSeq(inputFiles$normalFiles[i], formatName)
	}
	for(i in 1:nTumorRun) {
		tumorList[[i]] = readSeq(inputFiles$tumorFiles[i], formatName)
	}
	if(nNormalRun > 1) {
		normalSeq = CombineReadsAcrossRuns(normalList)
	}
	else {
		normalSeq = normalList[[1]]
	}
	if(nTumorRun > 1) {
		tumorSeq = CombineReadsAcrossRuns(tumorList)
	}
	else {
		tumorSeq = tumorList[[1]]
	}
	return(list(normalSeq=normalSeq, tumorSeq=tumorSeq))
}

