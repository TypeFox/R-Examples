

############################################################


sumpreserving.rounding <- 
function (data, digits = 0 , preserve = TRUE ) 
{
    ism <- ( is.matrix(data) ) | ( is.data.frame(data) )
    if (ism) {
        DD <- ncol(data)
        err.r <- data.r <- matrix(0, nrow = nrow(data), ncol = DD)
        data.r[, 1] <- round(data[, 1], digits)
			err.r[, 1] <- data[, 1] - data.r[, 1]
			for (dd in 2:DD) {
				data.r[, dd] <- round(data[, dd] + rowSums(err.r[, 1:(dd - 
					1),drop=FALSE]), digits)
				err.r[, dd] <- data[, dd] - data.r[, dd]
					}
		if (is.data.frame(data)){
			data.r <- data.frame(data.r)
			colnames(data.r) <- colnames(data)
				}
	   if ( ! preserve ){ data.r <- round( data , digits ) }
			}
		else {
			DD <- length(data)
			err.r <- data.r <- rep(0, DD)
			data.r[1] <- round(data[1], digits)
				err.r[1] <- data[1] - data.r[1]
				for (dd in 2:DD) {
					data.r[dd] <- round(data[dd] + sum(err.r[1:(dd - 
						1)]), digits)
					err.r[dd] <- data[dd] - data.r[dd]
					}
	   if ( ! preserve ){ data.r <- round( data , digits ) }
    }
    return(data.r)
}