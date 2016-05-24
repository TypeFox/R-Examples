Sliding.Window <-
function(VLF, seqlength, n = 30){
	samples <- nrow(VLF)
	window <- matrix(0, nrow = samples, ncol = seqlength-n)
	colors <- c("blue", "red", "green", "purple")
	
	for(r in 1:samples){
		for(i in 1:(seqlength-n)){
			for(z in 0:(n-1)){
				window[r,i] <- window[r,i] + VLF[r,i+z]
			}
			window[r,i] = window[r,i]/n
		}
	}
	
	theAxis = (c(1:(seqlength-n))/(seqlength-n))*100
	
	plot(x = theAxis, window[1,], type = "l", ylab = "VLFs/Position", xlab = "Percentile Barcode Segment", main = bquote("Sliding Window Analysis of ntVLFs (Window, N =" ~ .(n)~")"), col = colors[1])
	t = 2
	while(t <= samples){
		lines(x = theAxis, window[t,], type = "l", col = colors[t])
		t <- t + 1
	}
	if(samples > 1){
		legend("topright", legend = rownames(VLF), col = colors[1:samples], lty = 1)
	}
}
