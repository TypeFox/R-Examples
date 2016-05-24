LoadBeatEDFPlus <- function(HRVData, RecordName, RecordPath = ".", annotationType ="QRS", verbose = NULL) {
#-------------------------------
# Loads beats from an EDF+ file
#-------------------------------
#	RecordName -> record containing values
#	RecordPath -> path
#	typeAnnotation -> type of the annotation to read
#-------------------------------

	dir = getwd()
	setwd(RecordPath)

	VerboseDebug = FALSE

	if (!is.null(verbose)) {
        	cat("  --- Warning: deprecated argument, using SetVerbose() instead ---\n    --- See help for more information!! ---\n")
        	SetVerbose(HRVData, verbose)
    	}

	if (HRVData$Verbose) {
    	cat("** Loading beats positions for record:", RecordName,"**\n")
    	cat("   Path:", RecordPath, "\n")
	}

	con = file(RecordName,"rb")

	#version
	version=c()
	for(i in 1:8)
	{
		value = readBin(con,"integer",n=1,size=1,signed=FALSE)
		version=c(version,value)
	}

	if (VerboseDebug) {
		cat("------------- version:",rawToChar(as.raw(version)),"\n")
	}

	#patient identification
	patient=c()
	for(i in 1:80)
	{
		value = readBin(con,"integer",n=1,size=1,signed=FALSE)
		patient=c(patient,value)
	}

	if (VerboseDebug) {
		cat("------------- patient:",rawToChar(as.raw(patient)),"\n")
	}

	#recording identification
	recording=c()
	for(i in 1:80)
	{
		value = readBin(con,"integer",n=1,size=1,signed=FALSE)
		recording=c(recording,value)
	}

	if (VerboseDebug) {
		cat("------------- recording:",rawToChar(as.raw(recording)),"\n")
	}

	#date
	startdate=c()
	for(i in 1:8)
	{
		value = readBin(con,"integer",n=1,size=1,signed=FALSE)
		startdate=c(startdate,value)
	}

	if (VerboseDebug) {
		cat("------------- startdate:",rawToChar(as.raw(startdate)),"\n")
	}

	#time
	starttime=c()
	for(i in 1:8)
	{
		value = readBin(con,"integer",n=1,size=1,signed=FALSE)
		starttime=c(starttime,value)
	}

	if (VerboseDebug) {
		cat("------------- starttime:",rawToChar(as.raw(starttime)),"\n")
	}

	dateAux = intToCharacter(startdate)
	dateAux = gsub("\\.","-",dateAux)

	time = intToCharacter(starttime)
	time = gsub("\\.",":",time)

	if (HRVData$Verbose) {	
		cat("   Date: ",dateAux, "\n")
		cat("   Time: ",time, "\n")
	}
	datetimeinfo = paste(dateAux,time,sep = " ")
	HRVData$datetime=datetimeinfo

	#number of bytes in header
	numberBytes=c()
	for(i in 1:8)
	{
		value = readBin(con,"integer",n=1,size=1,signed=FALSE)
		numberBytes=c(numberBytes,value)
	}

	if (VerboseDebug) {
		cat("------------- numberBytes:",rawToChar(as.raw(numberBytes)),"\n")
	}

	#reserved
	reserved=c()
	for(i in 1:44)
	{
		value = readBin(con,"integer",n=1,size=1,signed=FALSE)
		reserved=c(reserved,value)
	}

	#number of data records
	numberDataRecords=c()
	for(i in 1:8)
	{
		value = readBin(con,"integer",n=1,size=1,signed=FALSE)
		numberDataRecords=c(numberDataRecords,value)
	}
	numberDataRecords = as.integer(intToCharacter(numberDataRecords))

	if (VerboseDebug) {
		cat("------------- numberDataRecords:",numberDataRecords,"\n")
	}

	#duration of a data record in seconds
	durationDataRecord=c()
	for(i in 1:8)
	{
		value = readBin(con,"integer",n=1,size=1,signed=FALSE)
		durationDataRecord=c(durationDataRecord,value)
	}

	if (VerboseDebug) {
		cat("------------- durationDataRecord:",rawToChar(as.raw(durationDataRecord)),"\n")
	}

	#number of signals in data record
	signalsDataRecord=c()
	for(i in 1:4)
	{
		value = readBin(con,"integer",n=1,size=1,signed=FALSE)
		signalsDataRecord=c(signalsDataRecord,value)
	}
	numberSignals = as.integer(intToCharacter(signalsDataRecord))

	if (VerboseDebug) {
		cat("------------- numberSignals:",numberSignals,"\n")
	}

	#Signals
	#label
	labels=c()
	limit = numberSignals * 16
	for(i in 1:limit)
	{
		value = readBin(con,"integer",n=1,size=1,signed=FALSE)
		labels=c(labels,value)
	}
	signalLabelList = stringToStringList(intToCharacter(labels))

	if (VerboseDebug) {
		cat("------------- labels:",labels,"\n")
		cat("------------- signalLabelList:",signalLabelList,"\n")
	}

	if (HRVData$Verbose) {
		cat("   Labels: ")
		for (i in 1:length(signalLabelList)) {
			cat(signalLabelList[i])
			if (i==length(signalLabelList)) {
				cat("\n")
			} else {
				cat(", ")
			}
		}
	}

	for(i in 1:length(signalLabelList))
	{
		if(signalLabelList[i]=="EDF")
		{
			if(i<length(signalLabelList) && signalLabelList[i+1]=="Annotations")
			{
				signalLabelList[i+1]=paste(signalLabelList[i],signalLabelList[i+1],sep=" ")
				signalLabelList = signalLabelList[-i]
			}
			break
		}
	}

	#transducer types
	types=c()
	limit = numberSignals * 80
	for(i in 1:limit)
	{
		value = readBin(con,"integer",n=1,size=1,signed=FALSE)
		types=c(types,value)
	}

	#physical dimension
	dimension=c()
	limit = numberSignals * 8
	for(i in 1:limit)
	{
		value = readBin(con,"integer",n=1,size=1,signed=FALSE)
		dimension=c(dimension,value)
	}

	#physical minimum
	minimum=c()
	limit = numberSignals * 8
	for(i in 1:limit)
	{
		value = readBin(con,"integer",n=1,size=1,signed=FALSE)
		minimum=c(minimum,value)
	}

	#physical maximum
	maximum=c()
	limit = numberSignals * 8
	for(i in 1:limit)
	{
		value = readBin(con,"integer",n=1,size=1,signed=FALSE)
		maximum=c(maximum,value)
	}

	#digital minimum
	dminimum=c()
	limit = numberSignals * 8
	for(i in 1:limit)
	{
		value = readBin(con,"integer",n=1,size=1,signed=FALSE)
		dminimum=c(dminimum,value)
	}

	#digital maximum
	dmaximum=c()
	limit = numberSignals * 8
	for(i in 1:limit)
	{
		value = readBin(con,"integer",n=1,size=1,signed=FALSE)
		dmaximum=c(dmaximum,value)
	}

	#prefiltering
	prefiltering=c()
	limit = numberSignals * 80
	for(i in 1:limit)
	{
		value = readBin(con,"integer",n=1,size=1,signed=FALSE)
		prefiltering=c(prefiltering,value)
	}

	#number of samples in each data record
	samples=c()
	limit = numberSignals * 8
	for(i in 1:limit)
	{
		value = readBin(con,"integer",n=1,size=1,signed=FALSE)
		samples=c(samples,value)
	}
	numberSamplesList = stringToIntList(intToCharacter(samples))

	#reserved
	reserved2=c()
	limit = numberSignals * 32
	for(i in 1:limit)
	{
		value = readBin(con,"integer",n=1,size=1,signed=FALSE)
		reserved2=c(reserved2,value)
	}

	if (HRVData$Verbose) {
		cat("   Reading data: this operation may take a few seconds...\n")
	}

	#Read annotations
	beats=c()

	for(i in 1:numberDataRecords)
	{
		for(j in 1:numberSignals)
		{
			if( signalLabelList[j] != "EDF Annotations")
			{
				limit = numberSamplesList[j]
				for(k in 1:limit)
				{
					value = readBin(con,"integer",n=1,size=2,signed=FALSE)
				}
			}
			else
			{
				limit = numberSamplesList[j] *2
				value=1
				while(value!=0)
				{
					value = readBin(con,"integer",n=1,size=1,signed=FALSE)
					limit=limit-1
				}

				seconds=c()
				annotation=c()

				value = readBin(con,"integer",n=1,size=1,signed=FALSE)
				limit=limit-1
				while(value!=0)
				{
					if(value==21)
					{
						while(value!=20)
						{
							value = readBin(con,"integer",n=1,size=1,signed=FALSE)
							limit=limit-1
						}
					}
					if(value!=20)
					{
						seconds=c(seconds,value)
					}
					else
					{
						value = readBin(con,"integer",n=1,size=1,signed=FALSE)
						limit=limit-1
						while(value!=0)
						{
							if(value==20)
							{
								if(intToCharacter(annotation)==annotationType)
								{
									seconds=intToCharacter(seconds)
									beats=c(beats,as.numeric(substr(seconds,2,nchar(seconds))))
								}
								annotation=c()
							}
							else
							{
								annotation=c(annotation,value)
							}
							value = readBin(con,"integer",n=1,size=1,signed=FALSE)
							limit=limit-1
						}

						seconds=c()
						annotation=c()
					}
					value = readBin(con,"integer",n=1,size=1,signed=FALSE)		
					limit=limit-1
				}
				for(k in 1:limit)
				{
					value = readBin(con,"integer",n=1,size=1,signed=FALSE)
				}
			}
		}
	}

	close(con)
	setwd(dir)

	if (length(beats)==0) {
		stop(" --- ERROR: Beats not found in file ",RecordName," ---\n   --- Quitting now!! ---\n")
	}

	HRVData$Beat = data.frame(Time = beats)

	if (HRVData$Verbose) {
        	cat("   Number of beats:", length(beats), "\n")
    	}

	return(HRVData)
}
