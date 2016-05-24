CreateHRVData <-
function(Verbose=FALSE) {
# ----------------------
# Creates the data model
# ----------------------
#	Verbose -> TRUE for Verbose mode

   	if (Verbose) {
      	cat("** Creating data model **\n")
   	}

   	HRVData <- list()
   	HRVData$Beat=NULL
		# Beat annotations: dataframe containing:
		#	Time -> sequence of beats (seconds)
		#	niHR -> non interpolated instantaneous frequency (beats/min)
		#	RR ->  instantaneous period (mseg.)
		
	HRVData$Episodes=NULL	
		# Episodes annotations: dataframe containing:
		#	InitTime -> initial time for the episodes (seconds)
		# 	Type -> type of episodes
		#  	Duration -> duration of episodes
		#  	Value -> numerical value for episodes
		
	HRVData$Freq_HR=NULL # Heart Rate interpolation frequency (default: 4 Hz)

	HRVData$HR=NULL # Interpolated Heart Rate 
	
	HRVData$Ext="hrv" # Extension for files
	
	HRVData$datetime=NULL # Date and time for the data

   	HRVData$FreqAnalysis <- list() # Frequency analysis associated to the data

   	HRVData$TimeAnalysis<-list() # Time analysis associated to the data
    
	HRVData$NonLinearAnalysis<-list() # Non Linear analysis associated to the data

   	HRVData$Verbose = Verbose
	
	if (HRVData$Verbose) {
		cat("   Data model created\n")
	}
	return(HRVData)
}

