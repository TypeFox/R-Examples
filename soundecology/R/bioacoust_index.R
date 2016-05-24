#Bioacoustic index
#
#Inspired from the "bioacoustic index" from the paper:
# Boelman NT, Asner GP, Hart PJ, Martin RE. 2007. Multi-trophic invasion
#   resistance in Hawaii: bioacoustics, field surveys, and airborne
#   remote sensing. Ecol Applications 17(8):2137-44.
#
# Inspired on Matlab code provided by NT Boelman. 
# Boelman et al. 2007 used min_freq=2000, max_freq=8000, fft_w=512
# Several parts where changed, in particular log math, so this won't be
# directly comparable to the original code in the paper.
#
# Requires: tuneR, seewave

bioacoustic_index <- function(soundfile, min_freq = 2000, max_freq = 8000, fft_w = 512){
	#Get sampling rate
	samplingrate <- soundfile@samp.rate
	freq_per_row = 10
	wlen = samplingrate/freq_per_row
	
	#Get Nyquist frequency in Hz
	nyquist_freq <- samplingrate/2
	
	if (max_freq > nyquist_freq) {
		cat(paste("\n ERROR: The maximum acoustic frequency that this file can use is ", nyquist_freq, "Hz. But the script was set to measure up to ", max_freq, "Hz.\n\n", sep = ""))
		break
		}
	
	#Stereo file
	if (soundfile@stereo == TRUE) {
		cat("\n This is a stereo file. Results will be given for each channel.\n")
		left <- channel(soundfile, which = c("left"))
		right <- channel(soundfile, which = c("right"))
		
		#Get values
		cat("\n Calculating index. Please wait... \n\n")
		spec_left <- spectro(left, f = samplingrate, wl = fft_w, plot = FALSE, dB = "max0")$amp
		spec_right <- spectro(right, f = samplingrate, wl = fft_w, plot = FALSE, dB = "max0")$amp
		#Clear from memory
		rm(left, right)
		
		#Get average in time
		specA_left <- apply(spec_left, 1, meandB)
		specA_right <- apply(spec_right, 1, meandB)
		
		#How much Hz are covered per row
		rows_width = length(specA_left) / nyquist_freq
		
		min_row = min_freq * rows_width
		max_row = max_freq * rows_width
		
		#Select rows
		specA_left_segment <- specA_left[min_row:max_row]
		specA_right_segment <- specA_right[min_row:max_row]
		
		freq_range <- max_freq - min_freq
		freqs <- seq(from = min_freq, to = max_freq, length.out = length(specA_left_segment))
		
		specA_left_segment_normalized <- specA_left_segment - min(specA_left_segment)
		specA_right_segment_normalized <- specA_right_segment - min(specA_right_segment)
		
		#left_area <- trapz(freqs, specA_left_segment_normalized)
		left_area <- sum(specA_left_segment_normalized * rows_width)
		right_area <- sum(specA_right_segment_normalized * rows_width)
		
		#cat("\n")
		cat("  Bioacoustic Index:\n")
		
		cat("   Left channel: ")
		cat(left_area)
		cat("\n   Right channel: ")
		cat(right_area)
		cat("\n\n")
	} else 
	{
		cat("\n This is a mono file.\n")
		#Get left channel
		left<-channel(soundfile, which = c("left"))
		
		#Get values
		cat("\n Calculating index. Please wait... \n\n")
		spec_left <- spectro(left, f = samplingrate, wl = fft_w, plot = FALSE, dB = "max0")$amp
		#Clear from memory
		rm(left)
		
		#Get average in time
		specA_left <- apply(spec_left, 1, meandB)
		
		#How much Hz are covered per row
		rows_width = length(specA_left) / nyquist_freq
		
		min_row = min_freq * rows_width
		max_row = max_freq * rows_width
		
		#Select rows
		specA_left_segment <- specA_left[min_row:max_row]
		freq_range <- max_freq - min_freq
		freqs <- seq(from = min_freq, to = max_freq, length.out = length(specA_left_segment))
		
		specA_left_segment_normalized <- specA_left_segment - min(specA_left_segment)
		
		#left_area <- trapz(freqs, specA_left_segment_normalized)
		left_area <- sum(specA_left_segment_normalized * rows_width)
		
		cat("  Bioacoustic Index: ")
		cat(left_area)
		cat("\n\n")
		right_area <- NA
	}
	invisible(list(left_area = left_area, right_area = right_area))
}
