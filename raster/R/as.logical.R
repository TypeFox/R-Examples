# Author: Robert J. Hijmans
# Date: November 2009
# Version 0.9
# Licence GPL v3


setMethod('as.logical', signature(x='Raster'), 
function(x, filename='', ...) {
	
	if (nlayers(x) > 1) {
		out <- brick(x, values=FALSE)
	} else {
		out <- raster(x)
	}
	
	if (canProcessInMemory(x, 2)){
		
		x <- getValues(x)
		x[] <- as.logical(x)
		out <- setValues(out, x)
		if (filename != '') {
			out <- writeRaster(out, filename, datatype='INT2S', ...)
		}
		return(out)
		
	} else {
		if (filename == '') {
			filename <- rasterTmpFile()					
		}
		
		out <- writeStart(out, filename=filename, ...)
		tr <- blockSize(x)
		pb <- pbCreate(tr$n, ...)	
		for (i in 1:tr$n) {
			v <- as.logical ( getValuesBlock(x, row=tr$row[i], nrows=tr$nrows[i] ) )
			out <- writeValues(out, v, tr$row[i])
			pbStep(pb, i) 
		} 
		pbClose(pb)			
		out <- writeStop(out)		
		return(out)
	}
}
)
