LoadBeatWFDB <- function (HRVData, RecordName, RecordPath = ".", annotator = "qrs", verbose = NULL) {
#-------------------------------
# Loads beats from a WFDB file
#-------------------------------
#	RecordName -> record containing values
#	RecordPath -> path
#	annotator -> type of file
#-------------------------------

	samplingFrequency = ""

	if (!is.null(verbose)) {
        	cat("  --- Warning: deprecated argument, using SetVerbose() instead ---\n    --- See help for more information!! ---\n")
        	SetVerbose(HRVData, verbose)
    	}
    	if (HRVData$Verbose) {
        	cat("** Loading beats positions for record:", RecordName,"**\n")
    	}

	dir = getwd() 
	setwd(RecordPath)

	auxHeader = readLines(paste(RecordName,".hea",sep=""),1)
	splitAuxHeader = strsplit(auxHeader," ")

	if(length(splitAuxHeader[[1]])>2)
		samplingFrequency = splitAuxHeader[[1]][3]
	else
		samplingFrequency = "250"

	samplingFrequency = as.numeric(samplingFrequency)


	#Read binary file
	con = file(paste(RecordName,".",annotator,sep=""),"rb")

	counter=1
	acumulator=0
	beats=c()
	repeat
	{
		value = readBin(con,"integer",n=1,size=2,signed=FALSE)

		binaryValue = intToBinary(value)

		while(length(binaryValue)<16)
			binaryValue = c(0,binaryValue)

		code = NULL
		time = NULL

		for(i in 1:6)
		{
			code=c(code,binaryValue[i])
		}

		for(i in 7:16)
		{
			time=c(time,binaryValue[i])
		}

		code=binaryToInt(code)
		time=binaryToInt(time)

		if(code==0 && time==0)
			break
		else
		{
			if(code==1)
			{
				acumulator=acumulator+time
				timeInSeconds = acumulator/samplingFrequency
				beats=c(beats,timeInSeconds)
				counter=counter+1
			}
			else
			{
				if(code==63)
				{
					jump = time%/%2 + time%%2
					for(i in 1:jump)
						value = readBin(con,"integer",n=1,size=2,signed=FALSE)
				}
				else
				{
					if(code==59 && time==0)
					{
						for(i in 1:2)
							value = readBin(con,"integer",n=1,size=2,signed=FALSE)
					}
					else
						if(code!=60 && code!=61 && code!=62 && code!=22 && code!=0)
						{
							acumulator=acumulator+time
						}
				}
			}
		}
	}

	HRVData$Beat = data.frame(Time = beats)

       	HRVData = LoadHeaderWFDB(HRVData, RecordName, RecordPath=".")

	if (HRVData$Verbose) {
        	cat("   Number of beats:", length(beats), "\n")
    	}

	close(con)
	setwd(dir)
	return(HRVData)
}
