#Sound Raster
	# Function to create an ASCII raster file from the spectrogram data. This
#  ASCII file can be directly opened in ArcGIS and most other GIS software.

sound_raster <- function(wavfile = NA, wav_directory = NA, max_freq = 10000, no_cores = 1){

  if (is.na(wavfile) == TRUE && is.na(wav_directory) == TRUE){
    stop(" You have to provide a filename, in the argument wavfile, or a directory as the argument wav_directory.\n\n")
  }
  
  if (is.na(wavfile) == FALSE && is.na(wav_directory) == FALSE){
  	stop(" You have to provide a single argument for which files to use, either wavfile or wav_directory.\n\n")
  }

  max_freq <- as.numeric(max_freq)

  write_ascii <- function(wavfile = NA, max_freq = max_freq, inCluster = FALSE){
  	#If launched in cluster, require the package for each node created
  	if (inCluster == TRUE){
  		require(soundecology)
  		}
	  	
  		soundfile <- readWave(wavfile)
	  	
	  	#Get Nyquist frequency in Hz
	  	nyquist_freq <- (soundfile@samp.rate/2)
	  	
	  	if (max_freq > nyquist_freq) {
	  		stop(paste("\n WARNING: The maximum acoustic frequency that this file can use is ", nyquist_freq, "Hz. But the script was set to measure up to ", max_freq, "Hz.\n\n", sep=""))
	  	}
	  	
	  	#window length for the spectro and spec functions
	  	#to keep each row every 10Hz
	  	#Frequencies and seconds covered by each
	  	freq_per_row = 10
	  	wlen = soundfile@samp.rate/freq_per_row
	  	
	  	#Stereo file
	  	if (soundfile@stereo == TRUE) {
	  		cat("\n This is a stereo file. A raster will be written for each channel. Please wait...\n")
	  		left <- channel(soundfile, which = c("left"))
	  		right <- channel(soundfile, which = c("right"))
	  		rm(soundfile)
	  		
	  		#matrix of values
	  		specA_left <- spectro(left, f = soundfile@samp.rate, wl = wlen, plot = FALSE)$amp
	  		specA_right <- spectro(right, f = soundfile@samp.rate, wl = wlen, plot = FALSE)$amp
	  		
	  		rm(left, right)
	  		
	  		#left channel
	  		max_freq_row <- max_freq/10
	  		sound_map_left <- specA_left[1:max_freq_row,]
	  		
	  		nrows <- dim(sound_map_left)[1]
	  		ncols <- dim(sound_map_left)[2]
	  		
	  		#Create container array for db values
	  		sound_map <- array(dim = dim(sound_map_left))
	  		
	  		#Using k to write the new array with lower freqs in the bottom
	  		k=dim(sound_map)[1]
	  		for (i in seq(from = 1, to = dim(sound_map)[1], by = 1)){
	  			for (j in seq(from = 1, to = dim(sound_map)[2], by = 1)){
	  				sound_map[k,j] = round(sound_map_left[i,j], digits = 4)
	  			}
	  			k = k-1
	  		}
	  		
	  		#ascii filename
	  		ascii_raster_left <- paste(strsplit(wavfile,".wav"), "_wav_left.asc", sep = "")
	  		#Write header
	  		cat(paste("ncols         ", ncols, "\nnrows         ", nrows, "\nxllcorner     0.0\nyllcorner     0.0\ncellsize      1\nNODATA_value  -9999\n", sep = ""), file = ascii_raster_left, append = FALSE)
	  		
	  		#Write data
	  		write.table(sound_map, file = ascii_raster_left, append = TRUE, row.names = FALSE, col.names = FALSE, sep = " ")
	  		
	  		
	  		
	  		#Right channel
	  		sound_map_right <- specA_right[1:max_freq_row,]
	  		
	  		#Create container array for db values
	  		sound_map <- array(dim = dim(sound_map_right))
	  		
	  		#Using k to write the new array with lower freqs in the bottom
	  		k <- dim(sound_map)[1]
	  		for (i in seq(from = 1, to = dim(sound_map)[1], by = 1)){
	  			for (j in seq(from = 1, to = dim(sound_map)[2], by = 1)){
	  				sound_map[k,j] <- round(sound_map_right[i,j], digits = 4)
	  			}
	  			k <- k-1
	  		}
	  		
	  		#ascii filename
	  		ascii_raster_right <- paste(strsplit(wavfile,".wav"), "_wav_right.asc", sep = "")
	  		#Write header
	  		cat(paste("ncols         ", ncols, "\nnrows         ", nrows, "\nxllcorner     0.0\nyllcorner     0.0\ncellsize      1\nNODATA_value  -9999\n", sep=""), file = ascii_raster_right, append=FALSE)
	  		
	  		#Write data
	  		write.table(sound_map, file = ascii_raster_right, append = TRUE, row.names = FALSE, col.names = FALSE, sep = " ")
	  		
	  		cat(paste("\n  ASCII raster files: \n    ", ascii_raster_left, "\n    ", ascii_raster_right, "\n\n", sep = ""))
	  		
	  	} else {
	  		
	  		cat("\n This is a mono file. A single raster will be written for this file. Please wait...\n")
	  		left<-channel(soundfile, which = c("left"))
	  		rm(soundfile)
	  		
	  		#matrix of values
	  		specA_left <- spectro(left, f = soundfile@samp.rate, wl = wlen, plot = FALSE)$amp
	  		
	  		rm(left)
	  		
	  		#left channel
	  		max_freq_row <- max_freq/10
	  		sound_map_left <- specA_left[1:max_freq_row,]
	  		
	  		nrows <- dim(sound_map_left)[1]
	  		ncols <- dim(sound_map_left)[2]
	  		
	  		#Create container array for db values
	  		sound_map <- array(dim = dim(sound_map_left))
	  		
	  		#Using k to write the new array with lower freqs in the bottom
	  		k=dim(sound_map)[1]
	  		for (i in seq(from = 1, to = dim(sound_map)[1], by = 1)){
	  			for (j in seq(from = 1, to = dim(sound_map)[2], by = 1)){
	  				sound_map[k,j] = round(sound_map_left[i,j], digits = 4)
	  			}
	  			k = k-1
	  		}
	  		
	  		#ascii filename
	  		ascii_raster_left <- paste(strsplit(wavfile,".wav"), "_wav_left.asc", sep = "")
	  		#Write header
	  		cat(paste("ncols         ", ncols, "\nnrows         ", nrows, "\nxllcorner     0.0\nyllcorner     0.0\ncellsize      1\nNODATA_value  -9999\n", sep=""), file = ascii_raster_left, append = FALSE)
	  		
	  		#Write data
	  		write.table(sound_map, file = ascii_raster_left, append = TRUE, row.names = FALSE, col.names = FALSE, sep = " ")
	  		
	  		cat(paste("\n  ASCII raster file: \n    ", ascii_raster_left, "\n\n", sep = ""))
	  		
	  	}
  }

  
  if (is.na(wavfile) == FALSE){
  	write_ascii(wavfile = wavfile, max_freq = max_freq)
  }else{
  	
  	if (file.access(wav_directory) == -1) {
  		stop(paste("The directory specified does not exist or this user is not autorized to read it:\n    ", wav_directory))
  	}
  	
  	#How many cores this machine has?
  	#require(parallel)
  	thismachine_cores <- detectCores()
  	
  	if (no_cores == 0){
  		stop("Number of cores can not be 0.")
  	}else if (no_cores < -1){
  		stop("Number of cores can not be negative.")
  	}else if (no_cores == "max"){
  		no_cores = thismachine_cores
  	}else if (no_cores == -1){
  		no_cores = thismachine_cores - 1
  	}else if (no_cores > thismachine_cores){
  		#Don't try to use more than the number of cores in the machine
  		warning(paste(" The number of cores to use can not be more than the
					  cores in this computer: ", detectCores()), immediate. = TRUE)
  		no_cores <- thismachine_cores
  	}
  	
  	wav_files <- dir(path = wav_directory, pattern = "wav$", ignore.case = TRUE, full.names = TRUE)
  	
  	if (length(wav_files) == 0) {
  		stop(paste(" No .wav files were found in the directory ", wav_directory, sep = ""))
  	}
  	
  	cat(paste("There are ", length(wav_files), " wav files.\n Please wait...\n", sep = ""))
  	
  	#Use parallel?
  	if (no_cores > 1){
  		no_files <- length(wav_files)
  		if (no_cores > no_files){
  			no_cores <- no_files
  			cat("\n The number of cores to use has been reduced because there are less files than cores available\n")
  		}
  		
  		cat(paste("\n Running on ", no_files, " files using ", no_cores, " cores", "\n\n", sep = ""))
  		
  		#Start timer
  		time0 <- proc.time()
  		
  		#create parallel cluster
  		cl <- makeCluster(no_cores, type = "PSOCK")
  		res <- parLapply(cl, wav_files, write_ascii, max_freq = max_freq, inCluster = TRUE)
  		
  		#Stop timer
  		time1 <- proc.time() - time0
  		cat(paste("\nProcess took: ", round(time1["elapsed"], 2), " seconds.\n\n", sep = ""))
  		
  		#pause to allow all to end
  		Sys.sleep(1)
  		
  		#Stop cluster
  		stopCluster(cl)
  		
  	}else{
  		
  		#Start timer
  		time0 <- proc.time()
	  		
	  	for (soundfile in wav_files){
	  		cat(paste("\n File: ", soundfile, sep=""))
	  		write_ascii(wavfile = soundfile, max_freq = max_freq)
	  	}
  		
  		#Stop timer
  		time1 <- proc.time() - time0
  		cat(paste("\nProcess took: ", round(time1["elapsed"], 2), " seconds.\n\n", sep = ""))
  		
  	}
  }
}
