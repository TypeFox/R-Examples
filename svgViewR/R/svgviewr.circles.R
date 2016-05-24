svgviewr.circles <- function(x, file=NULL, y=NULL, col=NULL, col.fill="black", 
	col.stroke="black", z.index=0, layer="", label="", r=1, lwd=2, opacity.stroke=1, 
	opacity.fill=1, append=TRUE, tag.name="point"){

	return_val <- svgviewr.points(x, file, y, type="p", col, col.fill, 
		col.stroke, z.index, layer, label, cex=r, lwd, opacity.stroke, 
		opacity.fill, append, tag.name="circle")

	return_val
}