energy.int <-
function(lim, k, A, rho) {
### internal function for energy
	
	bin.edges <- seq(lim[1], lim[2], by=1)
	bins <- tail(bin.edges, -1)-0.5
	wb.dist <- dweibull(bins, shape=k, scale=A)
	return(sum(0.5 * 8760 * rho * sum(bins^3 * wb.dist) / 10^3))
}
