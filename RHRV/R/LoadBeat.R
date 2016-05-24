LoadBeat <- function(fileType, HRVData, Recordname, RecordPath = ".", annotator = "qrs", scale = 1, datetime = "1/1/1900 0:0:0", annotationType = "QRS", verbose = NULL) {
#-------------------------------
# Loads beats from a specific file
#-------------------------------
#	fileType -> type of the file
#	RecordName -> record containing values
#	RecordPath -> path
#-------------------------------

	toret=""
	if(fileType == "WFDB")
	{
		toret = LoadBeatWFDB(HRVData,Recordname,RecordPath,annotator,verbose)
	}
	else
	{
		if(fileType == "Ascii")
		{
			toret = LoadBeatAscii(HRVData, Recordname, RecordPath, scale, datetime, verbose)
		}
		else
		{
			if(fileType == "RR")
			{
				toret = LoadBeatRR(HRVData, Recordname, RecordPath, scale, datetime, verbose)
			}	
			else
		 	{
				if(fileType == "Polar")
				{
					toret = LoadBeatPolar(HRVData, Recordname, RecordPath, verbose)
				}
				else
				{
					if(fileType == "Suunto")
					{
						toret = LoadBeatSuunto(HRVData, Recordname, RecordPath, verbose)
					}
					else
					{
						if(fileType == "EDFPlus")
						{
							toret = LoadBeatEDFPlus(HRVData, Recordname, RecordPath, annotationType, verbose)
						}
						else
						{
							print("Error: unknown file type")
						}
					}
				}
			}
		}
	}
	return(toret)
}
