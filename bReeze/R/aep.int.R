aep.int <-
function(wb.par, lim, pc, rho.pc, op, rho, avail) {
### internal function for calculation of annual energy production
	
	class.edges <- seq(lim[1], lim[2], by=1)
	classes <- tail(class.edges, -1)-0.5
	wb.dist <- dweibull(classes, shape=wb.par$k, scale=wb.par$A)
	pc.int <- spline(x=pc$v, y=pc$P, method="natural", xout=classes)[[2]]
	pc.int[pc.int<0] <- 0
		
	aep <- sum(wb.dist * pc.int) * op * rho/rho.pc  * avail / 1000
	aep[is.na(aep)] <- 0
	
	return(aep)
}
