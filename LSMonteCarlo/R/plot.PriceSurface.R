plot.PriceSurface <-
function(x, color=divPalette(800, "RdBu"), ...){
	persp(as.numeric(rownames(x)), as.numeric(colnames(x)), x, phi=25, theta=-45, ltheta=120, lphi=0, col=color,  ticktype="detailed", ylab="Strike", xlab="Volatility", zlab="Price", nticks=8, cex.lab=1.4, cex.axis=0.9)
}
