di.enm <- function(rn, x, y){
  
  if (!inherits(x, "Spatial")){
    if (!missing(rn)) {
      x <- gsub("_", paste("_", rn, "_", sep = ""), x)
    }
    x <- read.asciigrid(x)
  }
  
  if (!inherits(y, "Spatial")){
    if (!missing(rn)) {
      y <- gsub("_", paste("_", rn, "_", sep = ""), y)
    }
    y <- read.asciigrid(y)
  }
  
  xx <- slot(x, "data")
  yy <- slot(y, "data")
  
	# standardize probability surfaces
	# ------------------------------
	xSUM <- sum(xx, na.rm = TRUE)
	xx <- xx / xSUM
	ySUM <- sum(yy, na.rm = TRUE)
	yy <- yy / ySUM
	
	# Schoeners D (Schoener, 1968; Warren, Glor & Turelli, 2008)
	# ----------------------------------------------------------
	D <- 1 - 0.5 * sum(abs(xx - yy), na.rm = TRUE)
	
	# Hellingers Distance:
	# ---------------------
	H <- sqrt(sum((sqrt(xx) - sqrt(yy))^2, na.rm = TRUE))
	# I <- 1 - 0.5 * H -> error in Warren, Glor and Turelli 
	# (2008, Evolution 62:2868-2883)
	I <- 1 - H^2 * 0.5 # <- corrected I
  
  # both statistics range betweeen 0 (no overlap) and 1 (niches are identical)

	c(D = D, I = I)	
}