#' @export 
precintcon.deciles.analysis <- function(...) {
	
	l <- list(...)
	
	if (length(l) > 0) {
		
		varl <- as.list(match.call()[1:length(l)+1])
		
		quantl <- mapply(function(x, n) {
		
			if (is.element("precintcon.daily", class(x))) {
				x <- as.precintcon.monthly(x)
			} else if (!is.element("precintcon.monthly", class(x))) {
				stop("invalid input data. It should be either of type precintcon.daily or precintcon.monthly.")
			}
		
			q <- quantile(x[,ncol(x)], seq(.1, 1, by=.1))
			q <- as.vector(q)
			
			data <- data.frame(dataset = paste(n, sep = ""), 
				D1 = q[1], D2 = q[2], D3 = q[3], D4 = q[4], D5 = q[5],
				D6 = q[6], D7 = q[7], D8 = q[8], D9 = q[9], D10 = q[10])
			
			class(data) <- c("data.frame", "precintcon.deciles")
		
			return(data)
		}, l, varl, MoreArgs = list(), SIMPLIFY = FALSE)

		quantl <- do.call(rbind.data.frame, quantl)
		
		return(quantl)
		
	} else {
		stop("empty input data in precintcon.deciles.analysis function.")
	}
}