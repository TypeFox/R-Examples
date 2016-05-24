PlotSlope <- function(object, namemod = "default", 
	namex = "default", namey = "default", limitx = "default", 
	limity = "default") {
	yini<-yend<-NULL
	pmatr <- object$Points
	nomY <- object$nomY
	nomX <- object$nomX
	X_1L <- object$X_1L
	X_1H <- object$X_1H
	
	if (object$orde == 2) {
		nam <- dimnames(object$simple_slope)[1]
		nam <- nam[[1]]
		r1 <- nam[1]
		r2 <- nam[2]
		
		xini <- rep(X_1L, 4)
		xend <- rep(X_1H, 4)
		fact <- c(5, 6)
		mat <- cbind(fact, xini, pmatr[, 1], xend, pmatr[, 2])
		mat <- as.data.frame(mat)
		names(mat) <- c("fact", "xini", "yini", "xend", "yend")
		p <- ggplot(mat, aes(x = xini, y = yini))
		p1 <- p + geom_segment(aes(xend = xend, yend = yend))
		p1 <- p1 + scale_x_continuous(nomX) + scale_y_continuous(nomY)
		p1 <- p1 + geom_point(size = 3, aes(shape = factor(fact))) + 
			geom_point(aes(x = xend, y = yend, shape = factor(fact)), 
				size = 3)
		
		if (length(namemod) == 1) {
			p1 <- p1 + scale_shape(name = "Moderator", breaks = c(5, 
				6), labels = c(r1, r2))
		}
		if (length(namemod) > 1) {
			if (length(namemod) != 2) {
				stop("length of namemod vector must be = 2")
			}
			p1 <- p1 + scale_shape(name = "Moderator", breaks = c(5, 
				6), labels = namemod)
		}
		
		if (namex != "default") {
			if (length(limitx) == 2) {
				p1 <- p1 + scale_x_continuous(namex, limits = limitx)
			}
			else {
				p1 <- p1 + scale_x_continuous(namex)
			}
			
		}
		
		if (namey != "default") {
			if (length(limity) == 2) {
				p1 <- p1 + scale_y_continuous(namey, limits = limity)
			}
			else {
				p1 <- p1 + scale_y_continuous(namey)
			}
		}
		
		
		
		return(p1)
	}
	
	if (object$orde == 3) {
		
		nam <- dimnames(object$simple_slope)[1]
		nam <- nam[[1]]
		r1 <- nam[1]
		r2 <- nam[2]
		r3 <- nam[3]
		r4 <- nam[4]
		
		
		xini <- rep(X_1L, 4)
		xend <- rep(X_1H, 4)
		fact <- c(5, 6, 7, 8)
		mat <- cbind(fact, xini, pmatr[, 1], xend, pmatr[, 2])
		mat <- as.data.frame(mat)
		names(mat) <- c("fact", "xini", "yini", "xend", "yend")
		p <- ggplot(mat, aes(x = xini, y = yini))
		p1 <- p + geom_segment(aes(xend = xend, yend = yend))
		p1 <- p1 + scale_x_continuous(nomX) + scale_y_continuous(nomY)
		p1 <- p1 + geom_point(size = 3, aes(shape = factor(fact))) + 
			geom_point(aes(x = xend, y = yend, shape = factor(fact)), 
				size = 3)
		if (length(namemod) == 1) {
			p1 <- p1 + scale_shape(name = "Moderators Combination", 
				breaks = c(5, 6, 7, 8), labels = c(r1, r2, r3, 
				r4))
		}
		if (length(namemod) > 1) {
			if (length(namemod) != 4) {
				stop("length of namemod vector must be = 4")
			}
			p1 <- p1 + scale_shape(name = "Moderators Combination", 
				breaks = c(5, 6, 7, 8), labels = namemod)
		}
		p2 <- p1
		
		if (namex != "default") {
			if (length(limitx) == 2) {
				p2 <- p2 + scale_x_continuous(namex, limits = limitx)
			}
			else {
				p2 <- p2 + scale_x_continuous(namex)
			}
			
		}
		
		if (namey != "default") {
			if (length(limity) == 2) {
				p2 <- p2 + scale_y_continuous(namey, limits = limity)
			}
			else {
				p2 <- p2 + scale_y_continuous(namey)
			}
		}
		
		return(p2)
	}
	
}