#Measure Signals
# Function to extract range and other statistics of a
# song, sillable, or call of a sound file.

measure_signals <- function(wavfile, wl = 512, min_freq = NA, max_freq = NA, min_time = NA, max_time = NA, plot_range = 50, dBFS_range = 30, sample_size = 1, resultfile = NA, channel = "left"){

	if (is.na(resultfile) == TRUE){
		stop(" You have to provide a filename to save the results in the argument resultfile.")
		}
	
	soundfile1 <- readWave(wavfile)

	dir_to_save <- paste(strtrim(wavfile, nchar(wavfile) - 4), sep = "")
	soundfile <- channel(soundfile1, which = channel)
	rm(soundfile1)
	
	if (is.na(min_time) == TRUE){
		min_time = 0
		}
	if (is.na(max_time) == TRUE){
		max_time = length(soundfile@left) / soundfile@samp.rate
		}
	if (is.na(min_freq) == TRUE){
		min_fre = 0
		}
	if (is.na(max_freq) == TRUE){
		max_freq = soundfile@samp.rate / 2
		}
	
	#create results file, add column names
	fileheader <- c("WAVFILE,CHANNEL,WL,SAMPLE,PEAK_DBFS,PEAK_TIME,PEAK_FREQ,SIGNAL_TIME_MIN,SIGNAL_TIME_MAX,SIGNAL_FREQ_MIN,SIGNAL_FREQ_MAX,SELECTION_BOX_MINX,SELECTION_BOX_MAXX,SELECTION_BOX_MINY,SELECTION_BOX_MAXY")
	cat(fileheader, file = resultfile, append = FALSE)
	
	if (file.access(dir_to_save) == 0) {
		unlink(dir_to_save, recursive = TRUE)
		}
	dir.create(dir_to_save)
	
	specdata <- spectro(soundfile, wl = wl, flim=c(min_freq, max_freq), tlim=c(min_time, max_time), collevels = seq(-(abs(plot_range)), 0, length.out=16))
	
	for (i in 1:sample_size){
		
		cat(paste("\nSample ", i, ":\n", sep=""))
		cat(" Click on the four corners of a rectangle that contains the signal of interest:\n")
		plotpoints <- rbind(locator(1, type="o"))
		plotpoints <- rbind(plotpoints, rbind(locator(1, type="o")))
		plotpoints <- rbind(plotpoints, rbind(locator(1, type="o")))
		plotpoints <- rbind(plotpoints, rbind(locator(1, type="o")))
		
		minx <- min(unlist(plotpoints[, 1]))
		maxx <- max(unlist(plotpoints[, 1]))
		miny <- min(unlist(plotpoints[, 2]))
		maxy <- max(unlist(plotpoints[, 2]))
		
		#Get index of rows and cols for the selected area
		#from http://stackoverflow.com/a/10160443
		minx_index <- sapply(minx, function(x)which.min(abs(x - specdata$time)))
		maxx_index <- sapply(maxx, function(x)which.min(abs(x - specdata$time)))
		miny_index <- sapply(miny, function(x)which.min(abs(x - specdata$freq)))
		maxy_index <- sapply(maxy, function(x)which.min(abs(x - specdata$freq)))
		
		#Subset the spectrogram
		selection_data <- specdata$amp[miny_index:maxy_index, minx_index:maxx_index]
		selection_time <- specdata$time[minx_index:maxx_index]
		selection_freq <- specdata$freq[miny_index:maxy_index]
		
		#Find the peak inside the selection
		maxpeak_index <- which(selection_data == max(selection_data), arr.ind = TRUE)
		maxpeak_time <- selection_time[maxpeak_index[2]]
		maxpeak_freq <- selection_freq[maxpeak_index[1]]
		
		
		dfreqdata <- dfreq (soundfile, tlim = c(minx, maxx), bandpass = c(miny * 1000, maxy * 1000), ylim = c(miny, maxy), plot = FALSE)
		
		dfreqdata_x <- dfreqdata[,1]
		dfreqdata_y <- dfreqdata[,2]
		
		dfreqdata_rel <- cbind(dfreqdata_x + minx, dfreqdata_y)
		
		#plot dfreq on spectrogram
		lines(dfreqdata_rel, col="black")
		
		#write dfreq results to file
		dfreq_file <- paste(dir_to_save, "/", basename(wavfile), ".", i, ".csv", sep="")
		cat('"time (secs)", "frequency (kHz)"\n', file = dfreq_file)
		write.table(dfreqdata_rel, file = dfreq_file, append = FALSE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep=",")
		
		#boxplot(apply(selection_data, MARGIN=2, mean))
		#boxplot(apply(selection_data, MARGIN=1, mean))
		
		#selection_data_box <- boxplot(as.vector(selection_data))
		cutoff_value <- max(selection_data) - abs(dBFS_range)
		
		selection_data2 <- selection_data
		selection_data2[selection_data < cutoff_value] <- NA
		
		minsignalx <- selection_time[min(which( !is.na(selection_data2), arr.ind=TRUE)[, 2])]
		maxsignalx <- selection_time[max(which( !is.na(selection_data2), arr.ind=TRUE)[, 2])]
		minsignaly <- selection_freq[min(which( !is.na(selection_data2), arr.ind=TRUE)[, 1])]
		maxsignaly <- selection_freq[max(which( !is.na(selection_data2), arr.ind=TRUE)[, 1])]
		
		rect(xleft = minsignalx, ybottom = minsignaly, xright = maxsignalx, ytop = maxsignaly, border = "red", lwd = 3)
		
		this_res <- paste(basename(wavfile), channel, wl, i, max(selection_data), maxpeak_time, maxpeak_freq, minsignalx, maxsignalx, minsignaly, maxsignaly, minx, maxx, miny, maxy, sep=",")
		cat("\n", file = resultfile, append = TRUE)	
		cat(this_res, file = resultfile, append = TRUE)
	
		cat(paste("\nSaved the results to the file ", resultfile, "\n", sep=""))
		}
	
	}
