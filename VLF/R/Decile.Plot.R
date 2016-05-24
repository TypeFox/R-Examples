Decile.Plot <-
function(VLF, seqlength){
	samples <- nrow(VLF)
	decile <- matrix(0, nrow = samples, ncol = 10)
	colors <- c("blue", "red", "green", "purple")
	
	points <- seqlength/10
	each <- c()
	for(i in 1:10){
		each[i] <- round(points*i)
	}
	
	for(r in 1:samples){
		i <- 0
		for(n in 1:each[1]){
			decile[r,1] <- decile[r,1] + VLF[r, n]
		}
		for(i in 1:9){
			for(n in each[i]:each[i+1]){
				decile[r,i+1] <- decile[r,i+1] + VLF[r, n]
			}
		}
	}
	if(samples > 1){
		legend = rownames(VLF)}
	else{legend = NULL}
	
	b = barplot(decile, ylab = "Number of VLFs", xlab = "5' - Decile BARCODE Segment - 3'", main = "Distribution of VLFs Across Barcode", legend = rownames(VLF), col = colors[1:samples], beside = TRUE)
	b
	axis(1, at = b[seq(1,10*samples, samples)], labels = 1:10)
}
