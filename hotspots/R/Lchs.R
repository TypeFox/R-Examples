Lchs <- function(x, ...) {
if (!inherits(x, "hotspots")) 
	stop("use only with \"hotspots\" objects")
	x.s <- summary(x)
	x.Lc <- Lc(x.s$x)
	plot(x.Lc, ...)
	points(1-x.s$percent_phs*0.01,1-x.s$percent_phs_sum*0.01, cex = 2,pch = 16) }
