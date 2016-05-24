#' @noRd
#' @name precintcon.Z
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com} 
#' @aliases precintcon.Z 
#' @title precintcon.Z 
#' @description precintcon.Z 
#' @usage precintcon.Z(H) 
#' @param H H index
#' @return The Z index. 
#' @seealso 
#' \code{\link{precintcon.ci.analysis}}
#' @keywords precipitation Z index 
precintcon.Z <- function(H) {
	
	c0 <- 2.515517
	c1 <- 0.802853
	c2 <- 0.010328
	
	d1 <- 1.432788
	d2 <- 0.189269
	d3 <- 0.001308
	
	result <- c()
	
	for (prob in H) {
		t = 0
		minus = 0.0
		SPI = 0.0
		
		if (prob > 0.5) {
			minus = 1.0;
			prob = 1.0 - prob;
		} else {
			minus = -1.0;
		}
		
		if (prob < 0.0)
			return (0.0);

		if (prob == 0.0)
			SPI <- (9999.0 * minus)
		else {
			t = sqrt(log(1.0 / (prob * prob)));
		
			SPI <- (minus
						* (t
							- ((((c2 * t) + c1) * t) + c0)
							/ ((((((d3 * t) + d2) * t) + d1) * t) + 1.0)))
		}	
		result <- c(result, SPI)
	}
	
	return(result)
}