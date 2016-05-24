kzs.2d <- function(y, x, smooth, scale, k = 1, edges = TRUE, plot = TRUE)
{	
	s <- matrix(0, nrow = nrow(x) - 1, ncol = ncol(x))
	for (i in 1:2) {
		s[,i] <- diff(sort(x[,i]))
		if(scale[i] > min(s[,i][s[,i] > 0]))     
			stop("Each 'scale' must be less than or equal to the minimum of the difference of consecutive X values") 
		if(scale[i] <= 0)
			stop("Each 'scale' must be a positive real number")
		if(smooth[i] >= (max(x[,i]) - min(x[,i])))
			stop("Each 'smooth' must be much less than the difference of the max and min X") 
		if(smooth[i] <= 0)
			stop("Each 'smooth' must be a positive real number")
	}
	h <- smooth/2
	x1r <- x[,1]
	x2r <- x[,2]
	for (i in 1:k) {   
		x1i <- x[,1]
		x2i <- x[,2]   
		maxx <- apply(x, 2, max)                    
		minx <- apply(x, 2, min)
		yvals <- y
		xks <- vector("list", 2)
		for (i in 1:2) {
			xks[[i]] <- seq(minx[i] - h[i], maxx[i] + h[i], scale[i])
		}
		xk <- do.call("expand.grid", xks)   
		colnames(xk) <- c("x1", "x2")
		yk <- numeric(nrow(xk))
		for (j in 1:nrow(xk)) {
			w1 <- abs(x[,1] - xk$x1[j])
			w2 <- abs(x[,2] - xk$x2[j])
			Ik1 <- which(w1 <= h[1])
			Ik2 <- which(w2 <= h[2])
			Ik <- intersect(Ik1, Ik2)
			YIk <- yvals[Ik]
			yk[j] <- mean(YIk)     
		}
		df <- data.frame(cbind(xk, yk))             
		data <- na.omit(df)
		x <- as.matrix(data[,-3])
		y <- as.matrix(data$yk)
      }    
	if (edges == FALSE){
		between <- function (x,y) { x >= range(y)[1] & x <= range(y)[2] } 
		btw <- data[between(data[,1], x1r) & between(data[,2], x2r),]
		data <- na.omit(btw)
	}
	if (plot == TRUE){
		plot(wireframe(yk ~ x1 * x2, data, drape = TRUE, colorkey = TRUE, scales = list(arrows = FALSE)))    
	}		
	return(data)
}