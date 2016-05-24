"rntransform" <-
		function(formula,data,family=gaussian) {
	
	if ( is(try(formula,silent=TRUE),"try-error") ) { 
		if ( is(data,"gwaa.data") ) data1 <- phdata(data)
		else if ( is(data,"data.frame") ) data1 <- data
		else stop("'data' must have 'gwaa.data' or 'data.frame' class")
		formula <- data1[[as(match.call()[["formula"]],"character")]] 
	}
	
	var <- ztransform(formula,data,family)
	out <- rank(var) - 0.5
	out[is.na(var)] <- NA
	mP <- .5/max(out,na.rm=T)
	out <- out/(max(out,na.rm=T)+.5)
	out <- qnorm(out)
	out
}


