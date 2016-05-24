antimeridian.box <- function(object, tick.spacing=20){
	
	round(min(as.numeric(rownames(object))),2) -> lon.min
	round(max(as.numeric(rownames(object))),2) -> lon.max

	(lon.max-360) -> last.west

	lab.left <- rev(seq(180,lon.min,by=-tick.spacing))
	lab.right <- seq(-180, last.west, by=tick.spacing)[-1]
	lab <- c(lab.left,lab.right)
	n <- length(lab)
	
	lab2 <- lab
	lab2[lab2<0] <- lab2[lab2<0] + 360
	
	box()
	axis(2)
	axis(1, at=lab2, labels=lab)
}