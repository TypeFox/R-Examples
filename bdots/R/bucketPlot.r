bucket.plot <- function(buckets, ylim) {
	yellow <- rgb(255, 255, 0, alpha = 70, maxColorValue = 255)
	gray <- "gray44"
	if(!is.null(buckets)) {
		for(i in 1:nrow(buckets)) {
			rect(buckets[i,1], ylim[1], buckets[i,2], ylim[2], col = yellow, border = NA)
			lines(c(buckets[i,1], buckets[i,1]), ylim, col = gray, lwd = 1)
			lines(c(buckets[i,2], buckets[i,2]), ylim, col = gray, lwd = 1)
		}
	}
}