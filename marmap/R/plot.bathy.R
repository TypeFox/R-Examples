plot.bathy <- function(x, image=FALSE, bpal=NULL, land=FALSE, deepest.isobath, shallowest.isobath, step, n=20, lwd=1, lty=1, col="black", default.col="white", drawlabels = FALSE, xlab="Longitude", ylab="Latitude", asp=1, ...){
	
	x->mat # S3 compatibility
	
	if (!missing("deepest.isobath")) {
		length(deepest.isobath) -> n.levels
	} else {
		n.levels <- 1
	}
	
	unique(as.numeric(rownames(mat))) -> lon
	unique(as.numeric(colnames(mat))) -> lat
	  
	if (land == FALSE) mat[mat >0] <- 0

	if (image == FALSE){

		if(n.levels == 1){
			if (!missing("deepest.isobath") | !missing("shallowest.isobath") | !missing("step")) {
				seq(deepest.isobath, shallowest.isobath, by=step) -> level.param
			} else {
				level.param <- pretty(range(mat,na.rm=TRUE),n=n)
			}
			contour(lon,lat,mat, levels=level.param, 
					lwd=lwd, lty=lty, col=col, drawlabels = drawlabels, 
					xlab=xlab, ylab=ylab, asp=asp,...)
			}
		
		if(n.levels > 1){
			seq(deepest.isobath[1], shallowest.isobath[1], by=step[1]) -> level.param
			contour(lon,lat,mat,levels=level.param,
					lwd=lwd[1], lty=lty[1], col=col[1], drawlabels = drawlabels[1], 
					xlab=xlab, ylab=ylab, asp=asp, ...)
			for(i in 2:n.levels){
				seq(deepest.isobath[i], shallowest.isobath[i], by=step[i]) -> level.param
				contour(lon,lat,mat,levels=level.param,
						lwd=lwd[i], lty=lty[i], col=col[i], drawlabels = drawlabels[i], add=TRUE)
				}	
			# box(); axis(1); axis(2)
			}
		}

	if (image == TRUE) {

		if (is.null(bpal)) {
			colorRampPalette(c("#245372","#4871D9","#7D86A1","white")) -> ramp
			bpal <- ramp(100)
			}
			
		if (is.list(bpal))	bpal <- palette.bathy(mat, layers = bpal, land=land, default.col=default.col)

		if (n.levels == 1){
			if (!missing("deepest.isobath") | !missing("shallowest.isobath") | !missing("step")) {
				seq(deepest.isobath, shallowest.isobath, by=step) -> level.param
			} else {
				level.param <- pretty(range(mat,na.rm=TRUE),n=n)
			}
			image(lon,lat,mat, col=bpal, xlab=xlab, ylab=ylab, asp=asp, ...)
			contour(lon,lat,mat, levels=level.param, 
					lwd=lwd, lty=lty, col=col, drawlabels = drawlabels, add=TRUE)
			}
		
		if (n.levels > 1){
			image(lon,lat,mat, col=bpal, xlab=xlab, ylab=ylab, asp=asp, ...)
			seq(deepest.isobath[1], shallowest.isobath[1], by=step[1]) -> level.param
			contour(lon,lat,mat,levels=level.param,
					lwd=lwd[1], lty=lty[1], col=col[1], drawlabels = drawlabels[1], add=TRUE)
			for(i in 2:n.levels){
				seq(deepest.isobath[i], shallowest.isobath[i], by=step[i]) -> level.param
				contour(lon,lat,mat,levels=level.param,
						lwd=lwd[i], lty=lty[i], col=col[i], drawlabels = drawlabels[i], add=TRUE)
				}	
			}
		}
}