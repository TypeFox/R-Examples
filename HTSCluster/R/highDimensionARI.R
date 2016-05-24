
highDimensionARI <- function (x, y, splits = 2, verbose = FALSE) {

    	x <- as.vector(x)
    	y <- as.vector(y)
	len <- length(x)
	if(splits > floor(len/2)) {
		return("Error: Too many splits");
	}
	breaks <- round(seq(from = 0, to = len, length = splits + 1))

	a <- b <- c <- d <- 0
	for(i in 1:splits) {
		assign(paste("x", i, sep = ""), x[(breaks[i]+1):(breaks[i+1])])
		assign(paste("y", i, sep = ""), y[(breaks[i]+1):(breaks[i+1])])
		
		get.x <- get(paste("x", i, sep = ""))
		get.y <- get(paste("y", i, sep = ""))
		xx <- outer(get.x, get.x, "==")
		yy <- outer(get.y, get.y, "==")
		upper <- row(xx) < col(xx)
		xx <- xx[upper]
		yy <- yy[upper]
    		a <- a + sum(as.numeric(xx & yy))
    		b <- b + sum(as.numeric(xx & !yy))
    		c <- c + sum(as.numeric(!xx & yy))
    		d <- d + sum(as.numeric(!xx & !yy))

		if(verbose == TRUE) print(paste("Diag: ", i, sep = ""));
		rm(get.x)
		rm(get.y)
		rm(xx)
		rm(yy)
	}

	for(i in 1:(splits-1)) {
		for(j in (i+1):splits) {

			get.x1 <- get(paste("x", i, sep = ""))
			get.x2 <- get(paste("x", j, sep = ""))
			get.y1 <- get(paste("y", i, sep = ""))
			get.y2 <- get(paste("y", j, sep = ""))

			xx <- outer(get.x1, get.x2, "==")
			yy <- outer(get.y1, get.y2, "==")
    			a <- a + sum(as.numeric(xx & yy))
    			b <- b + sum(as.numeric(xx & !yy))
    			c <- c + sum(as.numeric(!xx & yy))
    			d <- d + sum(as.numeric(!xx & !yy))
			if(verbose == TRUE) {print(paste(i, "-", j))}
		
			rm(get.x1)
			rm(get.x2)
			rm(get.y1)
			rm(get.y2)
			rm(xx)
			rm(yy)
		}
		rm(list = paste("x", i, sep = ""))
		rm(list = paste("y", i, sep = ""))
	}
	
    	ni <- (b + a)
    	nj <- (c + a)
    	abcd <- a + b + c + d
    	q <- (ni * nj)/abcd
    	(a - q)/((ni + nj)/2 - q)
}
