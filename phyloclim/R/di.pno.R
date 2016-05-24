di.pno <- function(x, y){
  
	# normalize probability surfaces
	# ------------------------------
	xSUM <- sum(x, na.rm = TRUE)
	x <- x / xSUM
	ySUM <- sum(y, na.rm = TRUE)
	y <- y / ySUM
	
	# Schoeners D (Schoener, 1968; Warren, Glor & Turelli, 2008)
	# ----------------------------------------------------------
	D <- 1 - 0.5 * sum(abs(x - y), na.rm = TRUE)
	
	# Hellingers Distance:
	# ---------------------
	H <- sqrt(sum((sqrt(x) - sqrt(y))^2, na.rm = TRUE))
	# I <- 1 - 0.5 * H -> error in Warren, Glor and Turelli 
	# (2008, Evolution 62:2868-2883)
	I <- 1 - H^2 * 0.5 # <- corrected I
  
  # both statistics range betweeen 0 (no overlap) and 1 (niches are identical)

	c(D = D, I = I)	
}