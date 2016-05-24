add.Edges <-
function(grid, point.a, point.b){
	#points: c(long, lat)
	if (point.a[1] > point.b[1]){
		#point.a always on the left (lower x-value)
		tmp <- point.a
		point.a <- point.b
		point.b <- tmp
	}
	
	grid[point.a[1], point.a[2]] <- 1
	grid[point.b[1], point.b[2]] <- 1

	height <- abs(point.a[2] - point.b[2])
	width <- abs(point.a[1] - point.b[1])

	if ((height==0)&&(width==0)){
		#Case 1: point.a is equal point.b
		#--> no edge added
		return(grid)
	}
	
	if ((height==0)||(width==0)){
		#Case 2: straight lines (same x- or y-value)
		if (width == 0){
              	grid[point.a[1],point.a[2]:point.b[2]] <- rep(1, length(point.a[2]:point.b[2]))
       	}
        	if (height == 0){
             		grid[point.a[1]:point.b[1],point.a[2]] <- rep(1, length(point.a[1]:point.b[1]))
        	}
		return(grid)
	}

	if (height < width){
		quotient <- width/height
		step_width <- (point.b[1]-point.a[1])/width
           	step_height <- (point.b[2]-point.a[2])/height
		
		sum <- quotient/2
		sumvec <- round(sum, digits=3)
		jump <- c(round(sum))
		for (k in 0:(height-1)){
			sum <- sum + quotient
			sumvec <- c(sumvec,round(sum, digits=3))
			jump <- c(jump, round(sum))	
		}
		
		lastjump <- 0
		for (h in 0:height){
			nextjump <- ifelse(sumvec[h+1]-floor(sumvec[h+1]) == 0.5, floor(sumvec[h+1]),round(sumvec[h+1]))
			for (w in lastjump:nextjump){
				if (((point.a[1]+step_width*w)==point.b[1])&&((point.a[2]+step_height*h)==point.b[2])){
					break
				} else {
					grid[point.a[1]+step_width*w, point.a[2]+step_height*h] <- 1
				}
			}
			lastjump <-ifelse(sumvec[h+1]-floor(sumvec[h+1]) == 0.5,ceiling(sumvec[h+1]),round(sumvec[h+1]))
		}
	} 

	if (height == width){
		point.tmp <- point.a
		step_width <- (point.b[1]-point.a[1])/width
           	step_height <- (point.b[2]-point.a[2])/height
		while(point.tmp[1] != point.b[1]){
			point.tmp[1] <- point.tmp[1]+step_width
			point.tmp[2] <- point.tmp[2]+step_height
			grid[point.tmp[1], point.tmp[2]] <- 1
		} 
	}
	if (height > width){
		quotient <- height/width
		step_width <- (point.b[1]-point.a[1])/width
           	step_height <- (point.b[2]-point.a[2])/height
		
		sum <- quotient/2
		sumvec <- round(sum, digits=3)
		for (k in 0:(width-1)){
			sum <- sum + quotient
			sumvec <- c(sumvec,round(sum, digits=3))
		}
		
		lastjump <- 0
		for (w in 0:width){
			nextjump <- ifelse(sumvec[w+1]-floor(sumvec[w+1]) == 0.5, floor(sumvec[w+1]),round(sumvec[w+1]))
			for (h in lastjump:nextjump){
				if (((point.a[1]+step_width*w)==point.b[1])&&((point.a[2]+step_height*h)==point.b[2])){
					break
				} else {
					grid[point.a[1]+step_width*w, point.a[2]+step_height*h] <- 1
				}
			}
			lastjump <-ifelse(sumvec[w+1]-floor(sumvec[w+1]) == 0.5,ceiling(sumvec[w+1]),round(sumvec[w+1]))
		}
	}
	return(grid)
}