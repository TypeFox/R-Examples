TestData <- function(x, theta, sse, repeatable) {
	if (missing(repeatable)) 
		repeatable = F
	if (missing(sse)) 
		sse = 0.03
	if (missing(theta)) 
		theta <- c(-1, 1, 1, -0.24, 6, 0.2, 13)
	if (length(theta) == 3) 
		theta <- c(theta, 0, 0, 0, 0)
	if (length(theta) == 5) 
		theta <- c(theta, 0, 0)
	if (repeatable) 
		set.seed(143)

	Len <- length(x)
	Noise <- rnorm(Len, 0, sqrt(sse))
	Y <- theta[1] + theta[2] * exp(-x * theta[3]) + theta[4] * (x - theta[5]) * H(x, 10, theta[5]) + theta[6] * (x - 
		theta[7]) * H(x, 10, theta[7])

	tmp <- list(call = match.call(), time = x, thrs = (Y + Noise), resid = Noise, fit = Y, thet = theta, sse = sse, 
		val = var(Noise), data = paste("Test_data_", sse, sep = ""))
	class(tmp) = "dark"
	tmp
}