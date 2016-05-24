getMediansSubset <-
function(FILEPATH,NUMRUNSPERSAMPLE,MEASURES,RESULTFILENAME,ALTFILENAME,OUTPUTFILECOLSTART,OUTPUTFILECOLEND)
{
	# FUNCTION AS PART OF SPARTAN VERSION 2
	# THIS FUNCTION NO LONGER OUTPUTS A FILE, INSTEAD IT RETURNS THE LIST OF MEDIANS
	
	RESULTS<-NULL

	for(i in 1:NUMRUNSPERSAMPLE)
	{
		FILEADDRESS = paste(FILEPATH,toString(i),"/",RESULTFILENAME,sep="")
		#print(FILEADDRESS)
	
		if(file.exists(FILEADDRESS))
		{
			# CHECK WHAT KIND OF INPUT FILE IS BEING DEALT WITH
			if(substr(RESULTFILENAME,(nchar(RESULTFILENAME)+1)-3,nchar(RESULTFILENAME))=="csv")
			{
				# DEALING WITH A CSV FILE
				# import model result
				if(OUTPUTFILECOLSTART>1)
				{
					import<-read.csv(FILEADDRESS,colClasses=c(rep('NULL',OUTPUTFILECOLSTART-1),rep(NA,OUTPUTFILECOLEND-OUTPUTFILECOLSTART+1)),header=TRUE,check.names=FALSE)
				}else
				{
					import<-read.csv(FILEADDRESS,colClasses=c(rep(NA,OUTPUTFILECOLEND)),header=TRUE,check.names=FALSE)
				}

				MODELRESULT<-data.frame(import,check.names=FALSE)

			}else if(substr(RESULTFILENAME,(nchar(RESULTFILENAME)+1)-3,nchar(RESULTFILENAME))=="xml")
			{
				if(requireNamespace("XML",quietly=TRUE))
				{
					MODELRESULT<-XML::xmlToDataFrame(FILEADDRESS)
				}
				else
				{
					print("The getMediansSubset function requires the XML package to be installed")
				}
			}
			
			if(nrow(MODELRESULT)>0)
			{
				MEDIANSFORALLMEASURES<-NULL
	
				# NOW GET THE MEDIANS FOR EACH MEASURE
				for(q in 1:length(MEASURES))
				{
					# STORE JUST THE RESULTS FOR THAT MEASURE (sapply CAN THEN BE USED)
					MEASURE_RESULT<-as.matrix(MODELRESULT[MEASURES[q]])
					MEASUREMEDIAN<-median(as.numeric(MEASURE_RESULT))
					MEDIANSFORALLMEASURES<-cbind(MEDIANSFORALLMEASURES,MEASUREMEDIAN[[1]])
				}
				RESULTS<-rbind(RESULTS,MEDIANSFORALLMEASURES)
			}
			else		
			{
				if(!is.null(ALTFILENAME))
				{
					# USE THE ALTERNATIVE ON THIS OCCASION
					# ASSUMES IN SAME FORMAT AS ORIGINAL
					FILEADDRESS = paste(FILEPATH,toString(i),"/",ALTFILENAME,sep="")
					# import model result - again checking input type
					
					if(substr(ALTFILENAME,(nchar(ALTFILENAME)+1)-3,nchar(ALTFILENAME))=="csv")
					{
						# DEALING WITH A CSV FILE
						com <- paste("cut -d, -f",OUTPUTFILECOLSTART,"-",OUTPUTFILECOLEND," ",FILEADDRESS,sep="")
						import<-read.csv(pipe(com),header=TRUE,check.names=FALSE)
						MODELRESULT<-data.frame(import,check.names=FALSE)
					}
					else if(substr(ALTFILENAME,(nchar(ALTFILENAME)+1)-3,nchar(ALTFILENAME))=="xml")
					{
						if(requireNamespace("XML",quietly=TRUE))
						{
							MODELRESULT<-XML::xmlToDataFrame(FILEADDRESS)
						}
						else
						{
							print("The getMediansSubset function requires the XML package to be installed")
						}
					}
					
					if(nrow(MODELRESULT)>0)
					{
						MEDIANSFORALLMEASURES<-NULL
		
						# NOW GET THE MEDIANS FOR EACH MEASURE
						for(q in 1:length(MEASURES))
						{
							# STORE JUST THE RESULTS FOR THAT MEASURE (AS, INCASE OF XML, THIS WILL NEED TO BE
							# TRANSFORMED TO A MATRIX
							# IT WOULD BE NICE TO SPECIFY THE XML FIELD TYPES BUT THIS IS DIFFICULT TO DO WHEN INPUT
							# NEEDS TO BE GENERIC
							MEASURE_RESULT<-as.matrix(MODELRESULT[MEASURES[q]])
							MEASUREMEDIAN<-median(as.numeric(MEASURE_RESULT),na.rm = TRUE)
							# MEASUREMEDIAN<-sapply(MEASURE_RESULT,median)
							MEDIANSFORALLMEASURES<-cbind(MEDIANSFORALLMEASURES,MEASUREMEDIAN[[1]])
						}
						RESULTS<-rbind(RESULTS,MEDIANSFORALLMEASURES)
					}
				}
			}
		}
		else
		{
			print(paste("File ",FILEADDRESS," does not exist",sep=""))
		}
	}

	if(!is.null(RESULTS))
	{
		if(nrow(RESULTS)>0)
		{
			colnames(RESULTS)<-MEASURES
		}
	}

	# Now we return this set of medians
	return(RESULTS)
}

