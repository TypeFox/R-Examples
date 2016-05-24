col2alpha <- function(color, alpha = 0.5) {
	
	# Makes sure there is one alpha for each color
	if (length(alpha)!=length(color)) alpha <- rep(alpha,length(color))
	
	out <- numeric(length(color))
	for (i in 1:length(color)) {
		x <- col2rgb(color[i])[, 1]
		out[i] <- rgb(x[1], x[2], x[3], 255 * alpha[i], maxColorValue = 255)
	}
	out
}