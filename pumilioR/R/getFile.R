getFile <- function(result, SoundID = NA, credentials = NA){
	#Function to download the file associated with a SoundID
	if (is.na(SoundID)){
		stop(" SoundID cannot be empty.")
	}
	
	soundfilePath <- unlist(result[result$SoundID==SoundID,]$FilePath)
	localfile <- basename(soundfilePath)

	if (.Platform$OS.type == "windows") {
		
		if (!is.na(credentials)){
			#Fix for Windows systems, rarely is CURL installed and the default way
			# to download is very limited. This uses Internet Explorer functions
			#setInternet2(TRUE) 
			#^^^ build of package throws an error since it doesnt exist in Linux
			# Added to instructions			
			
			#Seems only way to give creds is by adding to address
			soundfilePath <- gsub("http://", paste("http://", credentials, "@", sep = ""), soundfilePath)
		}
		#	download.file(soundfilePath, destfile = localfile, mode = "wb")
		#}else{
			download.file(url = soundfilePath, destfile = localfile, mode = "wb")
		#}
	}else{
		#Mac and Linux		
		if (!is.na(credentials)){
			download.file(url = soundfilePath, destfile = localfile, extra = paste("-u",  credentials, sep=" "), mode="wb", method = "curl")
		}else{
			download.file(url = soundfilePath, destfile = localfile, mode = "wb", method = "curl")
		}
	  }
	
	#Return df of sound data
	invisible(localfile)
}