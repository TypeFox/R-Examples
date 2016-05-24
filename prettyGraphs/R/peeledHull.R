peeledHull <-
function(data_matrix,x_axis=1,y_axis=2,percentage=1,col="black",lwd=3,lty=1){
	nsim <- length(data_matrix[,x_axis])
	#I would prefer to use while(TRUE){}, rather than repeat.
	repeat{
		hpts <- chull(data_matrix[,c(x_axis,y_axis)])
		#if data_matrix[-hpts,] is only 1 row, npts is NULL.
		npts <- nrow(data_matrix[-hpts,c(x_axis,y_axis)])
		if((npts/nsim < percentage) || is.null(npts)){ 
			break 
		}
		data_matrix <- data_matrix[-hpts,]
	}
	hpts <- c(hpts,hpts[1])
	lines(data_matrix[hpts,c(x_axis, y_axis)], lwd=(lwd*1.5), lty=lty,col="black")
	lines(data_matrix[hpts,c(x_axis, y_axis)], lwd=lwd, lty=lty,col=col)
}
