# Copyright (C) 2011 Pierrick Bruneau, see README for full notice

RGBtoLab <- function(filename, filterWhite=FALSE, addCoords=TRUE) {
	# converts the ppm RGB file into Lab format	
	# only Lab data (no (x,y)) and filter white points
	
	# read file and build RGB matrix
	im <- pixmap::read.pnm(filename)
	red <- pixmap::getChannels(im, colors="red")
	green <- pixmap::getChannels(im, colors="green")
	blue <- pixmap::getChannels(im, colors="blue")
	
	# algo uses values in [0,1], which is what getChannels returns.
	
	n1 <- dim(red)[1]
	n2 <- dim(red)[2]
	
	cols <- matrix(nrow=0, ncol=3)
	coords <- matrix(nrow=0, ncol=2)
	for(i in 1:n1) {
		for(j in 1:n2) {
			if(filterWhite) {
				condition <- ((red[i,j] < 0.99) || (green[i,j] < 0.99) || (blue[i,j] < 0.99))
			} else {
				condition <- TRUE
			}		
			
			# filter white points
			if(condition) {
				if(red[i,j] > 0.04045) {
					red[i,j] <- ((red[i,j] + 0.055) / 1.055) ^ (2.4)
				} else {
					red[i,j] <- red[i,j] / 12.92
				}
				red[i,j] <- red[i,j] * 100
				if(green[i,j] > 0.04045) {
					green[i,j] <- ((green[i,j] + 0.055) / 1.055) ^ (2.4)
				} else {
					green[i,j] <- green[i,j] / 12.92
				}
				green[i,j] <- green[i,j] * 100
				if(blue[i,j] > 0.04045) {
					blue[i,j] <- ((blue[i,j] + 0.055) / 1.055) ^ (2.4)
				} else {
					blue[i,j] <- blue[i,j] / 12.92
				}
				blue[i,j] <- blue[i,j] * 100
			
			
				# with rbind => simpler
				cols <- rbind(cols, c(red[i,j], green[i,j], blue[i,j]))
				coords <- rbind(coords, c(j,i))
			}
		}
	}
	
	# transform matrix
	trs <- matrix(nrow=3, ncol=3)
	trs[1,] <- c(0.4124, 0.3576, 0.1805)
	trs[2,] <- c(0.2126, 0.7152, 0.0722)
	trs[3,] <- c(0.0193, 0.1192, 0.9505)	
	
	#trs <- trs / 0.17697
	
	# obtain xyz data
	cols <- cols %*% t(trs)
	
	
	#cols <- cols / const
	cols[,1] <- cols[,1] / 95.05
	cols[,2] <- cols[,2] / 100
	cols[,3] <- cols[,3] / 108.90
	
	thres <- 0.008856
	
	fct <- function(x, thres) { if(x > thres) { return(x^(1/3)) } else { return(7.787 * x + (16/116)) } }
	#cols <- apply(cols, 1:2, fct, thres=thres)
	for(i in 1:length(cols[,1])) {
		for(j in 1:3) {
			cols[i,j] <- fct(cols[i,j], thres)
		}
	}
	
	
	lab <- matrix(nrow=dim(cols)[1], ncol=3)
	lab[,1] <- 116 * cols[,2] - 16
	lab[,2] <- 500 * (cols[,1] - cols[,2])
	lab[,3] <- 200 * (cols[,2] - cols[,3])

	colmax <- max(coords[,1])
	coords[,1] <- coords[,1] / colmax
	colmax <- max(coords[,2])
	coords[,2] <- coords[,2] / colmax

	if(addCoords) {
		return(cbind(lab, coords))
	} else {
		return(lab)
	}
}







