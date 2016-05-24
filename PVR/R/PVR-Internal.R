.poss <- function(x){
	
	N <- 0
	for(i in 1:ncol(x)){
		
		N <- N + i
	}
	return(N)
}

.testPosition <- function(coords) {
	
	cx <- mean(coords[,1])
	cy <- mean(coords[,2])
	ix <- cy
	iy <- cx
	
	if(cx < ix & cy > iy){
		
		return(1)
	} else {
		
		return(-1)
	}
	
}


.Area <- function(coords){
	
	P <- round(coords, 7)
	T <- (round(data.frame(x = c(0, 1), y = c(0, 1)), 7))
	Pt <- (round(data.frame(x = 0, y = 0), 7))
	
	#first (sub)polygon first point 
	Pol <- (round(data.frame(x = P[1, 1], y = P[1, 2]), 7))
	
	result <- 0
	
	c1 <- 1
	for(i in 2:(nrow(P)-1)){
		
		if(!.Intersect(P[i:(i+1),], T)$intersect){
			
			Pol[i, ] <- P[i, ]
			Pol <- na.exclude(Pol)
		} else{
			
			Pol[i, ] <- P[i, ]
			#Pol <- na.exclude(Pol)
			
			if(P[i, 1] != P[i, 2]){
				
				Pt$x <- .Intersect(P[i:(i+1),], T)$point$x
				Pt$y <- .Intersect(P[i:(i+1),], T)$point$y
				
				Pol[(nrow(Pol) + 1), 1] <- Pt$x
				Pol[(nrow(Pol)), 2] <- Pt$y
				
				Pol <- na.exclude(Pol)
				result <- result + (.testPosition(Pol) * areapl(as.matrix(Pol)))
				
			} else {
				
				Pt <- P[i, ]
				
				Pol[(nrow(Pol) + 1), 1] <- Pt$x
				Pol[(nrow(Pol)), 2] <- Pt$y
				
				Pol <- na.exclude(Pol)
				
				result <- result + (.testPosition(Pol) * areapl(as.matrix(Pol)))
			}
			
			Pol <- data.frame(x = Pt[1, 1], y = Pt[1, 2])
			Pol <- na.exclude(Pol)
			
			c1 <- c1 + 1
			tmpPt <- Pt
		}	
	}
	return(result)
}

.Intersect <- function(l1, l2){
	
	res <- list(point = data.frame(x = NA, y = NA), intersect = FALSE)
	
	l1x1 <- l1[1 ,1]
	l1y1 <- l1[1 ,2]
	l1x2 <- l1[2 ,1]
	l1y2 <- l1[2 ,2]
	
	l2x1 <- l2[1 ,1]
	l2y1 <- l2[1 ,2]
	l2x2 <- l2[2 ,1]
	l2y2 <- l2[2 ,2]
	
	l1 <- data.frame(x = seq(l1x1,l1x2, length.out = 500), y = seq(l1y1,l1y2, length.out = 500))
	l2 <- data.frame(x = seq(l2x1,l2x2, length.out = 500), y = seq(l2y1,l2y2, length.out = 500))
	l1Eq <- lm(l1[,2] ~ l1[,1])
	l2Eq <- lm(l2[,2] ~ l2[,1])
	
	if(is.na(l2Eq$coefficients[2]) & is.na(l1Eq$coefficients[2])){
		
		if(l2[1, 1] == l1[1, 1]){
			
			res$intersect <- TRUE
			res$point$x <- l1[1, 1]
			res$point$y <- max(c(min(l1[,2]),min(l2[,2])))
			
			return(res)
		} else
			
			return(res)
	}
	
	if(is.na(l1Eq$coefficients[2])){
		
		x <- l1[1, 1]
		y <- l2Eq$coefficients[2] * l1[1, 1] + l2Eq$coefficients[1] 
		
		res$intersect <- TRUE
		res$point$x <- x
		res$point$y <- y
		
		return(res)
	} else
	
	if(is.na(l2Eq$coefficients[2])){
		
		x <- l2[1, 1]
		y <- l1Eq$coefficients[2] * l2[1, 1] + l1Eq$coefficients[1] 
		
		res$intersect <- TRUE
		res$point$x <- x
		res$point$y <- y
		
		return(res)
	} else
	
	if(!is.na(l1Eq$coefficients[2]) & !is.na(l2Eq$coefficients[2])){
		
		A <- matrix(c(-1,-1,round(l1Eq$coefficients[2], 7),round(l2Eq$coefficients[2], 7), round(l1Eq$coefficients[1], 7), round(l2Eq$coefficients[1], 7)),ncol = 3)
		x <- det(cbind(A[ ,3], A[ ,1]))/(det(A[ ,1:2]))
		y <- det(A[ ,2:3])/det(A[ ,1:2])
		
		res$intersect <- TRUE
		res$point$x <- round(x, 6)
		res$point$y <- round(y, 6)
		
		if(is.nan(res$point$x) & is.nan(res$point$y)){
			
			res$intersect <- TRUE
			res$point$x <- l1[1, 1]
			res$point$y <- l1[1, 2]
			
		}
		D_l1 <- sqrt(((l1[1, 1] - l1[nrow(l1), 1])^2) + ((l1[1, 2] - l1[nrow(l1), 2])^2))
		D_I1 <- sqrt(((l1[1, 1] - res$point[1, 1])^2) + ((l1[1, 2] - res$point[1, 2])^2))
		D_I2 <- sqrt(((res$point[1, 1] - l1[nrow(l1), 1])^2) + ((res$point[1, 2] - l1[nrow(l1), 2])^2))
		
		if( D_l1 < D_I1 | D_l1 < D_I2){
			
			res$intersect <- FALSE
			res$point$x <- NA
			res$point$y <- NA
			
			return(res)
		}
		return(res)
	}
}