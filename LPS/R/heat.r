# Color palette generation for heatmap (2 border + 1 center colors, symetric biases)
# Author : Sylvain Mareschal <maressyl@gmail.com>
heat.lin <- function(n, part) {
	shape <- 0:(n-1)
	return(shape)	
}
heat.exp <- function(n, part, base=1.015) {
	shape <- base^(0:(n-1))
	if(part == 2) shape <- max(shape) - rev(shape)
	return(shape)	
}
heat <- function(colors=c("#8888FF", "#000000", "#FF4444"), n=256, shapeFun=heat.exp, ...) {
	# Checks
	n <- as.integer(n)
	if(length(colors) != 3) stop("'colors' must contain 3 color names or hexadecimal representations")
	if(n %% 2L == 1L)       stop("'n' must be a multiple of 2")
	
	# Convert into RGB matrix
	input <- col2rgb(colors)
	
	# Interpolate
	output <- matrix(as.integer(NA), nrow=nrow(input), ncol=n, dimnames=list(rownames(input), NULL))
	for(channel in rownames(input)) {
		# Shapes
		shape1 <- shapeFun(n=n/2L, part=1L, ...)
		shape2 <- shapeFun(n=n/2L, part=2L, ...)
			
		# Color values
		output[channel,] <- c(
			shape1 / max(shape1) * (input[channel,2L] - input[channel,1L]) + input[channel,1L],
			shape2 / max(shape2) * (input[channel,3L] - input[channel,2L]) + input[channel,2L]
		)
	}
	
	# Matrix to hexadecimal representation
	output <- rgb(t(output), maxColorValue=255L)
	
	return(output)
}

