# Author: Jacob van Etten jacobvanetten@yahoo.com
# International Rice Research Institute
# Date :  January 2009
# Version 1.0
# Licence GPL v3

setMethod('raster', signature(x='TransitionLayer'), 
		function(x, reduceMethod="NZcolMeans") {
		rs <- as(x,"RasterLayer")
		dataVector <- vector(length=ncell(x))
		m <- as(x,"sparseMatrix")
		if(reduceMethod == "colSums") dataVector <- colSums(m)
		if(reduceMethod == "rowSums") dataVector <- rowSums(m)
		if(reduceMethod == "colMeans") dataVector <- colMeans(m)
		if(reduceMethod == "rowMeans") dataVector <- rowMeans(m)
		if(reduceMethod == "NZcolMeans" | reduceMethod == "NZrowMeans"){
			mL <- as(m,"lMatrix")
			if(reduceMethod == "NZrowMeans") dataVector <- rowSums(m)/rowSums(mL)
			if(reduceMethod == "NZcolMeans") dataVector <- colSums(m)/colSums(mL)
		}
		rs <- setValues(rs, dataVector) 
		return(rs)
	}
)