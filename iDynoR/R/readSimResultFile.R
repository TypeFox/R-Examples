readSimResultFile <-
function(resultFolder,resultFileType,timePoint)
{
	opSys<-.Platform$OS.type
	if(opSys=="windows")
	{
		# REPLACE THE SLASH IN THE FILEPATH IF PRESENT
		resultFolder<-gsub("\\\\","/",resultFolder)
		
	}

	#### READ IN THE XML FILE
	doc <- xmlTreeParse(paste(resultFolder,"/",resultFileType,"/",resultFileType,"(",timePoint,").xml",sep=""))
	#### GET THE DATA
	xmlResultData<-xmlRoot(doc)

	return(xmlResultData)
}
