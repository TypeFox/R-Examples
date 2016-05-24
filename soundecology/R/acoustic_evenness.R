#Acoustic Evenness Index from Villanueva-Rivera \emph{et al.} 2011. 
# The ADI is calculated by dividing the spectrogram into bins (default 10) and taking the proportion of the signals in each bin above a threshold (default -50 dBFS). The ADI is the result of the Gini index applied to these bins.

acoustic_evenness <- function(soundfile, max_freq = 10000, db_threshold = "-50", freq_step = 1000){
	
	db_threshold <- as.numeric(db_threshold)
		
	#function that gets the proportion of values over a db
	# value in a specific band of frequencies. 
	# Frequency is in Hz
	getscore<-function(spectrum, minf, maxf, db, freq_row){
		miny<-round((minf)/freq_row)
		maxy<-round((maxf)/freq_row)
		
		subA=spectrum[miny:maxy,]
		
		index1<-length(subA[subA>db])/length(subA)
		
		return(index1)
	}
	
	#Some general values
	#Get sampling rate
	samplingrate <- soundfile@samp.rate
	
	#Get Nyquist frequency in Hz
	nyquist_freq <- samplingrate/2
	
	#window length for the spectro and spec functions
	#to keep each row every 10Hz
	#Frequencies and seconds covered by each
	freq_per_row = 10
	wlen = samplingrate/freq_per_row
	
	#Stereo file
	if (soundfile@stereo == TRUE) {
		cat("\n This is a stereo file. Results will be given for each channel.\n")
		left<-channel(soundfile, which = c("left"))
		right<-channel(soundfile, which = c("right"))
		rm(soundfile)
		
		#matrix of values
		cat("\n Calculating index. Please wait... \n\n")
		specA_left <- spectro(left, f = samplingrate, wl = wlen, plot = FALSE)$amp
		specA_right <- spectro(right, f = samplingrate, wl = wlen, plot = FALSE)$amp
		
		rm(left,right)
		
		if (max_freq > nyquist_freq) {
			cat(paste("\n WARNING: The maximum acoustic frequency that this file can use is ", nyquist_freq, "Hz. But the script was set to measure up to ", max_freq, "Hz. The value of max_freq was changed to ", nyquist_freq, ".\n\n", sep=""))
			max_freq <- nyquist_freq
		}
		
		Freq<-seq(from = 0, to = max_freq - freq_step, by = freq_step)
		
		#LEFT CHANNEL
		Score <- rep(NA, length(Freq))
		
		for (j in 1:length(Freq)) {
			Score[j]=getscore(specA_left, Freq[j], (Freq[j]+freq_step), db_threshold, freq_per_row)
		}
		
		left_vals=Score
		
		#RIGHT CHANNEL
		Score <- rep(NA, length(Freq))
		
		for (j in 1:length(Freq)) {
			Score[j]=getscore(specA_right, Freq[j], (Freq[j]+freq_step), db_threshold, freq_per_row)
		}
		
		right_vals = Score
		
# 		cat(" ==============================================\n")
# 		cat(paste(" Results (with a dB threshold of ", db_threshold, ")\n\n", sep=""))
		
		left_bandvals_return <- rep(NA, length(Freq))
		right_bandvals_return <- rep(NA, length(Freq))
		left_bandrange_return <- rep(NA, length(Freq))
		right_bandrange_return <- rep(NA, length(Freq))
		
# 		cat(" Proportion over threshold for each frequency band (in csv format): \n\n")
# 		cat("Frequency range (Hz), left channel proportion, right channel proportion\n")
		for (j in seq(length(Freq), 1, by = -1)) {
# 			cat(paste(Freq[j], "-", (Freq[j]+freq_step), ",", round(left_vals[j],6), ",", round(right_vals[j],6), "\n", sep=""))
			left_bandvals_return[j] = round(left_vals[j], 6)
			right_bandvals_return[j] = round(right_vals[j], 6)
			left_bandrange_return[j] = paste(Freq[j], "-", (Freq[j] + freq_step), " Hz", sep = "")
			right_bandrange_return[j] = paste(Freq[j], "-", (Freq[j] + freq_step), " Hz", sep = "")
		}
		
# 		cat("\n Plot of proportions in each band: \n\n")
# 		cat("  Left channel\n")
# 		cat("   Freq. range (Hz) |--------------------|\n")
		
		#printed in inverse order to keep the low frequencies in the bottom, like in a spectrogram
		for (j in seq(length(Freq), 1, by = -1)) {
			this_row_name <- paste(Freq[j], "-", (Freq[j] + freq_step), "", sep = "")
			this_row_size <- nchar(this_row_name)
			this_row_space <- 17 - this_row_size
			
			this_row_spaces = ""
			
			for (f in seq(1,this_row_space,by = 1)) {
				this_row_spaces = paste(this_row_spaces, " ", sep = "")
			}
			
# 			cat(paste("   ", this_row_name, this_row_spaces, "|", sep=""))
# 			temp_val=round(left_vals[j],2)*20
# 			if (temp_val>0){
# 				for (i in 1:temp_val) {
# 					cat("*")
# 				}
# 			}
# 			cat("\n")
# 			rm(temp_val)
		}
		
# 		cat("\n  Right channel\n")
# 		cat("   Freq. range (Hz) |--------------------|\n")
		
		#printed in inverse order to keep the low frequencies in the bottom, like in a spectrogram
		for (j in seq(length(Freq), 1, by = -1)) {
			this_row_name <- paste(Freq[j], "-", (Freq[j] + freq_step), "", sep = "")
			this_row_size <- nchar(this_row_name)
			this_row_space <- 17 - this_row_size
			
			this_row_spaces = ""
			
			for (f in seq(1,this_row_space, by = 1)) {
				this_row_spaces = paste(this_row_spaces, " ", sep="")
			}
			
# 			cat(paste("   ", this_row_name, this_row_spaces, "|", sep=""))
			
# 			temp_val=round(right_vals[j],2)*20
# 			if (temp_val>0){
# 				for (i in 1:temp_val) {
# 					cat("*")
# 				}
# 			}
# 			cat("\n")
# 			rm(temp_val)
		}
		
		#cat("\n")
		cat("  Acoustic Evenness Index: \n")
		cat(paste("   Left channel: ", round(Gini(left_vals), 6), "\n", sep=""))
		cat(paste("   Right channel: ", round(Gini(right_vals), 6), "\n\n", sep=""))
		left_gini_return = round(Gini(left_vals), 6)
		right_gini_return = round(Gini(right_vals), 6)
		
	} else 
	{
		cat("\n This is a mono file.\n")
		
		#matrix of values
		cat("\n Calculating index. Please wait... \n\n")
		specA_left <- spectro(soundfile, f = samplingrate, wl = wlen, plot = FALSE)$amp
		
		rm(soundfile)
		
		if (max_freq>nyquist_freq) {
			cat(paste("\n WARNING: The maximum acoustic frequency that this file can use is ", nyquist_freq, "Hz. But the script was set to measure up to ", max_freq, "Hz. The value of max_freq was changed to ", nyquist_freq, ".\n\n", sep=""))
			max_freq <- nyquist_freq
		}
		
		Freq<-seq(from = 0, to = max_freq - freq_step, by = freq_step)
		
		#Score=seq(from=0, to=0, length=length(Freq))
		Score <- rep(NA, length(Freq))
		
		for (j in 1:length(Freq)) {
			Score[j] = getscore(specA_left, Freq[j], (Freq[j] + freq_step), db_threshold, freq_per_row)
		}
		
		left_vals = Score
		
		
# 		cat(" ==============================================\n")
# 		cat(paste(" Results (with a dB threshold of ", db_threshold, ")\n\n", sep=""))
# 		
# 		cat(" Proportion over threshold for each frequency band (in csv format): \n\n")
# 		cat("Frequency range (Hz), proportion\n")
		
		#printed in inverse order to keep the low frequencies in the bottom, like in a spectrogram
		left_bandvals_return <- rep(NA, length(Freq))
		right_bandvals_return <- rep(NA, length(Freq))
		left_bandrange_return <- rep(NA, length(Freq))
		right_bandrange_return <- rep(NA, length(Freq))
		for (j in seq(length(Freq), 1, by = -1)) {
# 			cat(paste(Freq[j], "-", (Freq[j]+freq_step), ",", round(left_vals[j],6), "\n", sep=""))
			left_bandvals_return[j] = round(left_vals[j], 6)
			left_bandrange_return[j] = paste(Freq[j], "-", (Freq[j] + freq_step), " Hz", sep = "")
		}
		
# 		cat("\n Plot of proportions in each band: \n\n")
# 		cat("   Freq. range (Hz) |--------------------|\n")
		
		#printed in inverse order to keep the low frequencies in the bottom, like in a spectrogram
		for (j in seq(length(Freq), 1, by = -1)) {
			this_row_name <- paste(Freq[j], "-", (Freq[j] + freq_step), "", sep = "")
			this_row_size <- nchar(this_row_name)
			this_row_space <- 17 - this_row_size
			
			this_row_spaces = ""
			
			for (f in seq(1, this_row_space, by = 1)) {
				this_row_spaces = paste(this_row_spaces, " ", sep = "")
			}
			
# 			cat(paste("   ", this_row_name, this_row_spaces, "|", sep=""))
# 			temp_val=round(left_vals[j],2)*20
# 			if (temp_val>0){
# 				for (i in 1:temp_val) {
# 					cat("*")
# 				}
# 			}
# 			cat("\n")
# 			rm(temp_val)
		}
		
		#cat("\n")
		cat("  Acoustic Evenness Index: ")
		cat(paste(round(Gini(left_vals), 6), "\n", sep = ""))
		left_gini_return = round(Gini(left_vals), 6)
		right_gini_return = NA
	}
	invisible(list(aei_left = left_gini_return, aei_right = right_gini_return))
}
