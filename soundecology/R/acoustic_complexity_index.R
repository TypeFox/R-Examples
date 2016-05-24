#
#Acoustic Complexity Index
#From: N. Pieretti, A. Farina, D. Morri. 2011. A new methodology to infer
# the singing activity of an avian community: The Acoustic Complexity Index (ACI).
# Ecological Indicators 11: 868-873.
#
#Tested with SoundscapeMeter 1.0.14.05.2012, courtesy of A. Farina
#
acoustic_complexity <- function(soundfile, min_freq = NA, max_freq = NA, j = 5, fft_w = 512){
	
  if (is.na(max_freq)){
    max_freq <- soundfile@samp.rate / 2
  }
  
  if (is.na(min_freq)){
    min_freq <- 0
  }
  
	#function that gets the difference of values
	get_d <- function(spectrum, freq_row, min_col, max_col){
		D = 0
		for (k in min_col:(max_col - 1)) {
			D = D + abs(spectrum[freq_row,k] - spectrum[freq_row,k + 1])
			}
			 		
		return(D)
		}
	
	#Some general values
	#Get sampling rate
	samplingrate <- soundfile@samp.rate
	duration <- length(soundfile@left)/soundfile@samp.rate
	
	#Get Nyquist frequency in Hz
	nyquist_freq <- (samplingrate/2)
	if (max_freq>nyquist_freq) {
		cat(paste("\n WARNING: The maximum acoustic frequency that this file can use is ", nyquist_freq, "Hz. But the script was set to measure up to ", max_freq, "Hz. The value of max_freq was changed to ", nyquist_freq, ".\n\n", sep = ""))
		max_freq <- nyquist_freq
		#break
	}
	
	#window length for the spectro and spec functions
	wlen = fft_w
	
	#Stereo file
	if (soundfile@stereo == TRUE) {
		cat("\n This is a stereo file. Results will be given for each channel.\n")
		left <- channel(soundfile, which = c("left"))
		right <- channel(soundfile, which = c("right"))
				
		#matrix of values
		cat("\n Calculating index. Please wait... \n\n")
		spec_left <- spectro(left, f = samplingrate, wl = wlen, plot = FALSE, norm = TRUE, dB = NULL, scale = FALSE, wn = "hamming")
		
		specA_left <- spec_left$amp
		
		min_freq1k = min_freq/1000 
		max_freq1k = max_freq/1000 
		
		which_min_freq <- which(abs(spec_left$freq - min_freq1k)==min(abs(spec_left$freq - min_freq1k)))
		which_max_freq <- which(abs(spec_left$freq - max_freq1k)==min(abs(spec_left$freq - max_freq1k))) 
		  
		if (which_min_freq <1){
		  which_min_freq = 1
		}
		
		if (which_max_freq > dim(specA_left)[1]){
		  which_max_freq = dim(specA_left)[1]-1
		}
		
# 		cat(which_min_freq)
# 		cat(",")
# 		cat(which_max_freq)
# 		cat(",")
# 		cat(dim(specA_left))
		specA_left <- spec_left$amp[which_min_freq:which_max_freq,]
		rm(spec_left)
		
		spec_right <- spectro(right, f = samplingrate, wl = wlen, plot = FALSE, norm = TRUE, dB = NULL, scale = FALSE, wn = "hamming")
		
		specA_right <- spec_right$amp[which_min_freq:which_max_freq,]
		
		rm(spec_right)
		
		rm(left,right)
		
# 		specA_rows <- dim(specA_left)[1]
# 		specA_cols <- dim(specA_left)[2]
# 		
# 		freq_per_row <- specA_rows/nyquist_freq
# 				
# 		max_row <- round(max_freq * freq_per_row)
# 		
# 		specA_left <- specA_left[1:max_row,]
# 		specA_right <- specA_right[1:max_row,]
		specA_rows <- dim(specA_left)[1]
		specA_cols <- dim(specA_left)[2]
		
		fl <- rep(NA, specA_rows)
		delta_fl <- ( max_freq - min_freq ) / specA_rows
		delta_tk <- (length(soundfile@left)/soundfile@samp.rate) / specA_cols
				
		#m <- floor(duration / j)
		#q <- specA_rows
		no_j <- floor(duration / j)
		
		#Number of values, in each row, for each j period (no. of columns)
		I_per_j <- floor(j/delta_tk)
		
		ACI_left_vals <- rep(NA, no_j)
		ACI_fl_left_vector <- rep(NA, no_j)
		ACI_left_matrix <- data.frame(matrix(NA, nrow = specA_rows, ncol = no_j))
		
		ACI_right_vals <- rep(NA, no_j)
		ACI_fl_right_vector <- rep(NA, no_j)
		ACI_right_matrix <- data.frame(matrix(NA, nrow = specA_rows, ncol = no_j))
		
		#Left channel
		#For each frequency bin fl
		for (q_index in 1:specA_rows) {
			
			#For each j period of time
			for (j_index in 1:no_j) {
				min_col <- j_index * I_per_j - I_per_j + 1
				max_col <- j_index * I_per_j
								
				D <- get_d(specA_left, q_index, min_col, max_col)
				sum_I <- sum(specA_left[q_index,min_col:max_col])
				ACI_left_vals[j_index] <- D / sum_I
				ACI_left_matrix[q_index, j_index] <- D / sum_I
			}
			
			ACI_fl_left_vector[q_index] <- sum(ACI_left_vals)
			} 
		
		ACI_tot_left <- sum(ACI_fl_left_vector)
		
		#Right channel
		#For each frequency bin fl
		for (q_index in 1:specA_rows) {
			
			#For each j period of time
			for (j_index in 1:no_j) {
				min_col <- j_index * I_per_j - I_per_j + 1
				max_col <- j_index * I_per_j				
				
				D <- get_d(specA_right, q_index, min_col, max_col)
				sum_I <- sum(specA_right[q_index, min_col:max_col])
				ACI_right_vals[j_index] <- D / sum_I
				ACI_right_matrix[q_index, j_index] <- D / sum_I
			}
			
			ACI_fl_right_vector[q_index] <- sum(ACI_right_vals)
			} 
		
		ACI_tot_right <- sum(ACI_fl_right_vector)
		
		ACI_tot_left_by_min <- round((ACI_tot_left/duration) * 60, 2)
		ACI_tot_right_by_min <- round((ACI_tot_right/duration) * 60, 2)
		
		cat(paste("  Acoustic Complexity Index (total):\n", "   Left channel: ", sep=""))
		cat(ACI_tot_left)
		cat(paste("\n", "   Right channel: ", sep=""))
		cat(ACI_tot_right)
		cat("\n\n")
		if (duration > 60){
  		cat(paste("  Acoustic Complexity Index (by minute):\n", "   Left channel: ", sep=""))
  		cat(ACI_tot_left_by_min)
  		cat(paste("\n", "   Right channel: ", sep=""))
  		cat(ACI_tot_right_by_min)
  		cat("\n\n")
  		}
		
	} else 
	{
		cat("\n This is a mono file.\n")
		
		left<-channel(soundfile, which = c("left"))
				
		#matrix of values
		cat("\n Calculating index. Please wait... \n\n")
		spec_left <- spectro(left, f = samplingrate, wl = wlen, plot = FALSE, norm = TRUE, dB = NULL, scale = FALSE, wn = "hamming")
		
		specA_left <- spec_left$amp
		
		min_freq1k = min_freq/1000 
		max_freq1k = max_freq/1000 
		
		which_min_freq <- which(abs(spec_left$freq - min_freq1k)==min(abs(spec_left$freq - min_freq1k)))
		which_max_freq <- which(abs(spec_left$freq - max_freq1k)==min(abs(spec_left$freq - max_freq1k))) 
		
		specA_left <- specA_left[which_min_freq:which_max_freq,]
		rm(spec_left)
		
		rm(left)
		
		#LEFT CHANNEL
		specA_rows <- dim(specA_left)[1]
		specA_cols <- dim(specA_left)[2]
# 		
# 		freq_per_row <- specA_rows/nyquist_freq
# 				
# 		max_row <- round(max_freq * freq_per_row)
# 		
# 		specA_left <- specA_left[1:max_row,]
# 		specA_rows <- dim(specA_left)[1]
		
		fl <- rep(NA, specA_rows)
		delta_fl <- ( max_freq - min_freq ) / specA_rows
		delta_tk <- (length(soundfile@left)/soundfile@samp.rate) / specA_cols
		
		no_j <- floor(duration / j)
		#q <- specA_rows
		#m <- floor(duration / j)
		
		#Number of values, in each row, for each j period (no. of columns)
		I_per_j <- floor(j/delta_tk)
		
		ACI_left_vals <- rep(NA, no_j)
		ACI_fl_left_vector <- rep(NA, no_j)
		ACI_left_matrix <- data.frame(matrix(NA, nrow = specA_rows, ncol = no_j))
		
		ACI_right_vals <- rep(NA, no_j)
		ACI_fl_right_vector <- rep(NA, no_j)
		ACI_right_matrix <- data.frame(matrix(NA, nrow = specA_rows, ncol = no_j))
		
		#Left channel
		#For each frequency bin fl
		for (q_index in 1:specA_rows) {
			
			#For each j period of time
			for (j_index in 1:no_j) {
				min_col <- j_index * I_per_j - I_per_j + 1
				max_col <- j_index * I_per_j
								
				D <- get_d(specA_left, q_index, min_col, max_col)
				sum_I <- sum(specA_left[q_index, min_col:max_col])
				ACI_left_vals[j_index] <- D / sum_I
				ACI_left_matrix[q_index, j_index] <- D / sum_I
			}
			
			ACI_fl_left_vector[q_index] <- sum(ACI_left_vals)
		} 
		
		ACI_tot_left <- sum(ACI_fl_left_vector)
		ACI_tot_left_by_min <- round((ACI_tot_left/duration) * 60, 2)
		
		ACI_tot_right <- NA
		ACI_tot_right_by_min <- NA
		
		cat("  Acoustic Complexity Index (total): ")
		cat(ACI_tot_left)
		cat("\n\n")
		if (duration > 60){
  		cat("  Acoustic Complexity Index (by minute): ")
  		cat(ACI_tot_left_by_min)
  		cat("\n\n")
		  }
	}
	
	invisible(list(AciTotAll_left = ACI_tot_left, AciTotAll_right = ACI_tot_right, 
	               AciTotAll_left_bymin = ACI_tot_left_by_min, AciTotAll_right_bymin = ACI_tot_right_by_min,
				   #AciIfTotAll_left=ACIif_tot_left, AciIfTotAll_right=ACIif_tot_right, 
				   aci_fl_left_vals = ACI_fl_left_vector, aci_fl_right_vals = ACI_fl_right_vector,
				   #aci_if_left_vals=ACI_if_left_vector, aci_if_right_vals=ACI_if_right_vector,
				   aci_left_matrix = ACI_left_matrix, aci_right_matrix = ACI_right_matrix))
}
