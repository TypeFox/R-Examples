c_and_d <- function(s, t){
	c1.vec <- c2.vec <- c3.vec <- NULL
	d1.vec <- d2.vec <- d3.vec <- NULL
	for(c3 in 0:min(c(t, s))){
		for(d3 in 0:min(c(t - c3, s - c3))){
			for(c1 in 0:floor((t - c3 - d3)/2)){
				for(c2 in 0:floor((s - c3 - d3)/2)){
					if((t - 2*c1 - c3 - d3) %% 2 == 0 & (s - 2*c2 - c3 - d3) %% 2 == 0){
						d1 <- (t - 2*c1 - c3 - d3)/2
						d2 <- (s - 2*c2 - c3 - d3)/2
						c1.vec <- c(c1.vec, c1)
						c2.vec <- c(c2.vec, c2)
						c3.vec <- c(c3.vec, c3)
						d1.vec <- c(d1.vec, d1)
						d2.vec <- c(d2.vec, d2)
						d3.vec <- c(d3.vec, d3)
					}
				}
			}
		}
	}
	return(data.frame(c1 = c1.vec, c2 = c2.vec, c3 = c3.vec,
			d1 = d1.vec, d2 = d2.vec, d3 = d3.vec))
}