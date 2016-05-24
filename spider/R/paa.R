paa <- 
function(data, sppVector){
	singPoly <- function(vec){
		aa <- unique(vec)
		bb <- length(aa)
		cc <- ifelse(bb == 1, aa, "poly")
		cc
		}
	apply(data, MARGIN = 2, FUN = function(x) tapply(x, sppVector, singPoly))
}