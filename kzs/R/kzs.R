kzs <- function(y, x, smooth, scale, k = 1, edges = TRUE, plot = TRUE)
{	
	x <- sort(x) 
	dx <- diff(x)
	if(scale > min(dx[dx > 0]))
		stop("'scale' must be less than or equal to the minimum difference of consecutive X values")  
	if(smooth >= (max(x) - min(x)))
		stop("'smooth' must be much less than the difference of the max and min X values") 
	h <- smooth/2
	xrange <- range(x)
	for (i in 1:k) {
		xi <- as.vector(x)                      
		maxx <- max(xi)	                    
		minx <- min(xi)			        
		yvals <- y          
		xk <- seq(minx - h, maxx + h, scale)           
		yk <- numeric(length(xk))
       	for(j in 1:length(xk)) {
			w <- abs(xi - xk[j])     	         
			Ik <- which(w <= h)	
			size <- length(Ik)	
			Yik <- yvals[Ik]			
			yk[j] <- (1 / size) * sum(Yik)	  
		}                    
		df <- data.frame(cbind(xk, yk))             
		data <- na.omit(df)      
		x <- data$xk
		y <- data$yk            
	}              	
	if(edges == FALSE){
		edgs <- data[data$xk >= min(xrange) & data$xk <= max(xrange), ]   
		data <- as.data.frame(edgs)
	}
	if (plot == TRUE) {
		plot(data$yk ~ data$xk, xlab = "xk", ylab = "yk", xlim = c(min(xrange), max(xrange)), 
		     type = 'l', pch = 19, col = 'black')  
	}
	return(data)
}