# Lab and RGB functions : used to derive a percepually uniform
# color channel (instead of gray scales)
#
# This transform is based on ITU-R BT 709, using the D65
# white point reference. 
# See http://fr.wikipedia.org/wiki/Rec._709 for details.
# If RGB vals scale in [0,1], L should scale in [0,100], and (a,b) in [-110,110].
# also see http://www.easyrgb.com/index.php?X=MATH
RGB2Lab <- function(mat) {
	# cast mat in case of a single numeric vector
	mat <- matrix(mat, ncol=3)
	# input should be a matrix with R,G and B as columns in [0,1], and columnwise pixels as rows.
	# RGB -> XYZ
	thres1 <- 0.04045
	M <- c(0.412453, 0.357580, 0.180423,
		   0.212671, 0.715160, 0.072169,
		   0.019334, 0.119193, 0.950227)
	M <- matrix(M, nrow=3, byrow=TRUE)
	matthres <- mat > thres1
	mat <- matthres * ((mat + 0.055) / 1.055) ^ 2.4 + (!matthres) * mat / 12.92
	xyz <- mat %*% t(M)
	
	# XYZ -> Lab
	thres2 <- 0.008856
	xyz <- sweep(xyz, 2, c(0.950456,1,1.088754), "/")
	#yalone <- xyz[,2]
	#y3 <- yalone^(1/3)
	xyzthres <- xyz > thres2
	xyz <- xyzthres * xyz^(1/3) + (!xyzthres) * (7.787*xyz+16/116)
	#L <- xyzthres[,2] * (116*y3-16) + (!xyzthres[,2]) * (903.3*yalone)
	L <- 116 * xyz[,2] - 16
	a <- 500 * (xyz[,1] - xyz[,2])
	b <- 200 * (xyz[,2] - xyz[,3])
	return(cbind(L,a,b))
}

Lab2RGB <- function(mat) {
	# cast mat in case of a single numeric vector
	mat <- matrix(mat, ncol=3)
	# same spec as RGBtoLab
	# Lab -> XYZ
	thres1 <- 0.008856
	thres2 <- 0.0031308
	
	y <- (mat[,1] + 16) / 116
	x <- mat[,2] / 500 + y
	z <- y - mat[,3] / 200
	xyz <- cbind(x,y,z)
	xyzthres <- xyz > thres1
	xyz <- xyzthres * (xyz^3) + (!xyzthres) * ((xyz - 16 / 116) / 7.787)
	xyz <- sweep(xyz, 2, c(0.950456,1,1.088754), "*")
	
	# XYZ -> RGB
	Minv = c(3.240479,-1.537150,-0.498535,
            -0.969256, 1.875992, 0.041556,
             0.055648,-0.204043, 1.057311)
    Minv <- matrix(Minv, nrow=3, byrow=TRUE)
    RGB <- xyz %*% t(Minv)
	RGBthres <- RGB > thres2
	# manage NaN risk with x^a and a<1
	RGB[RGBthres] <- 1.055 * RGB[RGBthres]^(1/2.4) - 0.055
	RGB[!RGBthres] <- 12.92 * RGB[!RGBthres]
	# bound in case of very small absolute <0 values, and resp. for >1 values.
	RGB[RGB<0] <- 0
	RGB[RGB>1] <- 1
    return(RGB)
}


# utility function to probe a single color, specified as rgb or lab vector (no alpha supported) with values in [0,1]
drawSinglePatch <- function(vec, rgb=TRUE) {
	# if rgb is FALSE, lab is assumed.
	if (!rgb) vec <- Lab2RGB(vec)
	pushViewport(viewport())
	vp1 <- viewport(x=unit(0.5, "native"), y=unit(0.5, "native"), width=unit(0.4, "native"), height=unit(0.4, "native"))
	grid.rect(x=unit(0, "snpc"), y=unit(0, "snpc"), width=unit(1, "snpc"), height=unit(1, "snpc"), gp=gpar(col=rgb(vec[1],vec[2],vec[3]), fill=rgb(vec[1],vec[2],vec[3])), just=c("left", "bottom"), vp=vp1)	
	return(NULL)
}

# utility function to assess a pseudocolor sequence (interpolation between 2 RGB colors performed in Lab space)
drawGradient <- function(col1, col2, rgb=TRUE, numcols=200) {
	# generate a set of colors with linear interpolation in lab space
	if (rgb) {
		col1 <- RGB2Lab(col1)
		col2 <- RGB2Lab(col2)
	}
	
	Linter <- seq(from=col1[1], to=col2[1], length.out=numcols)
	ainter <- seq(from=col1[2], to=col2[2], length.out=numcols)
	binter <- seq(from=col1[3], to=col2[3], length.out=numcols)
	graddata <- cbind(Linter, ainter, binter)
	graddata <- Lab2RGB(graddata)
	coords <- seq(from=0, to=1, length.out=numcols)
	pushViewport(viewport())
	vp1 <- viewport(x=unit(0.5, "native"), y=unit(0.5, "native"), width=unit(0.8, "native"), height=unit(0.3, "native"))	
	grid.rect(x=unit(coords, "npc"), y=unit(rep(0, numcols), "npc"), width=unit(1/numcols, "npc"), height=unit(1, "npc"), gp=gpar(col=rgb(graddata[,1],graddata[,2],graddata[,3]), fill=rgb(graddata[,1],graddata[,2],graddata[,3])), just=c("left","bottom"), vp=vp1)
	return(NULL)
}



genHatchData <- function(rgbcol, rgbhatch, side=20) {
	rgbcol <- matrix(rgbcol, ncol=3)
	rgbhatch <- matrix(rgbhatch, ncol=3)
	step <- (side / 10) * 2
	hatchInds <- seq(from=(1+step), to=side, by=step)
	rgbpatch <- list()
	class(rgbpatch) <- "rgbpatch" # class rgbpatch can manage one to n image patches.
	
	datafill <- function(x) {
		mat <- matrix(x[1], nrow=side, ncol=side)
		mat[hatchInds,] <- mat[,hatchInds] <- x[2]
		return(mat)
	}

	rgbpatch$rchannel <- t(apply(cbind(rgbcol[,1], rgbhatch[,1]), 1, datafill))
	rgbpatch$gchannel <- t(apply(cbind(rgbcol[,2], rgbhatch[,2]), 1, datafill))
	rgbpatch$bchannel <- t(apply(cbind(rgbcol[,3], rgbhatch[,3]), 1, datafill))
	return(rgbpatch)
}




drawPatches <- function(data, patches, patchSize=0.05, alpha=0.5, patchNpix=20, highlight=numeric(0), labels=rep("1", length(highlight))) {
	if (class(patches) != "rgbpatch") stop("patches should be an rgbpatch object (see genHatchDatafor for reference)")
	nrows <- dim(data)[1]
	if (nrows != dim(patches$rchannel)[1]) stop("patches and data should reference the same number of elements")
	
	xlim <- c(min(data[,1]), max(data[,1]))
	ylim <- c(min(data[,2]), max(data[,2]))
	margin <- c(xlim[2]-xlim[1], ylim[2]-ylim[1]) * 0.05
	pushViewport(viewport(xscale=c(xlim[1]-margin[1],xlim[2]+margin[1]), yscale=c(ylim[1]-margin[2],ylim[2]+margin[2])))
	patchSize <- unit(patchSize, "snpc") # convert to normalized units
	if (length(highlight) == 0) {
		alphavec <- rep(alpha, nrows)
		labs <- rep("", nrows)
	} else {
		alphavec <- rep(alpha/20, nrows)
		alphavec[highlight] <- alpha
		labs <- rep("", nrows)
		labs[highlight] <- labels
	}
	alpha <- alphavec
	
	for (i in 1:nrows) {
		rowi <- rep(1:patchNpix, times=patchNpix)
		coli <- rep(1:patchNpix, each=patchNpix)
		curvp <- viewport(x=unit(data[i,1], "native"), y=unit(data[i,2], "native"), width=patchSize, height=patchSize, just=c("center", "center"))
		grid.rect(x=unit((coli-1)/patchNpix, "snpc"), y=unit((rowi-1)/patchNpix, "snpc"), width=unit(1/patchNpix, "snpc"), height=unit(1/patchNpix, "snpc"), gp=gpar(col=rgb(patches$rchannel[i,], patches$gchannel[i,], patches$bchannel[i,],alpha[i]), fill=rgb(patches$rchannel[i,], patches$gchannel[i,], patches$bchannel[i,], alpha[i])), just=c("left", "bottom"), vp=curvp)
		grid.text(labs[i], y=unit(-0.2, "snpc"), vp=curvp)
	}
	return(NULL)
}





animatePatches <- function(nsecs=5) {
	reftime <- getTimestamp()
	laptime <- 5 # seconds for a complete revolution
	initangle <- 0
	pushViewport(viewport())
	curvp <- viewport(x=unit(0.5 + 0.3*sin(initangle), "npc"), y=unit(0.5 + 0.3*cos(initangle), "npc"), width=unit(0.1, "npc"), height=unit(0.1, "npc"))
	grid.rect(name="myrect", x=unit(0.5, "npc"), y=unit(0.5, "npc"), width=unit(1.0, "npc"), height=unit(1.0, "npc"), gp=gpar(col="black", fill="black"), vp=curvp)
	while(getElapsed(reftime) < nsecs) {
		newangle <- initangle + (getElapsed(reftime) %% laptime) * (2*pi) / laptime
		curvp <- viewport(x=unit(0.5 + 0.3*sin(newangle), "npc"), y=unit(0.5 + 0.3*cos(newangle), "npc"), width=unit(0.1, "npc"), height=unit(0.1, "npc"))
		grid.edit("myrect", vp=curvp)
		Sys.sleep(1/30)
	}
	return(NULL)
}


plotGreyPatches <- function(data, maingrey, auxgrey=NULL, patchSize=0.05, alpha=0.5, highlight=numeric(0), labels=rep("1", length(highlight))) {
	# draw black-bordered grey patches. Grey levels in maingrey indicate some visual attribute.
	# auxgrey may indicate an additional attribute. Was originally designed for displaying projective distorsions (see semisupKernelPCA package).
	# better ways of displaying this may exist (see aupetit 2007), but didnt want to rely on color cues.
	# note that using grey levels as here is only for a rough estimation : due to contrast effects, accurate value and comparison judgements should not be made (see Ware) 
	# data frames and matrices work the same with apply : column-wise. just ensure correct conversions with character data frames.
	
	nrows <- dim(data)[1]
	d <- dim(data)[2]
	if (is.null(auxgrey)) {
		method <- "squarepatch"
	} else {
		method <- "trianglepatch"
	}
	
	xlim <- c(min(data[,1]), max(data[,1]))
	ylim <- c(min(data[,2]), max(data[,2]))
	margin <- c(xlim[2]-xlim[1], ylim[2]-ylim[1]) * patchSize
	pushViewport(viewport(xscale=c(xlim[1]-margin[1],xlim[2]+margin[1]), yscale=c(ylim[1]-margin[2],ylim[2]+margin[2])))
	margin <- unit(margin, "native")
	if (length(highlight) == 0) {
		alphavec <- rep(alpha, nrows)
		labs <- rep("", nrows)
	} else {
		alphavec <- rep(alpha/7, nrows)
		alphavec[highlight] <- alpha
		labs <- rep("", nrows)
		labs[highlight] <- labels
	}
	
	alpha <- alphavec
	
	squarepatch <- function(frame) {
		coords <- as.numeric(frame[1:d])
		grey <- as.numeric(frame[d+1])
		alpha <- as.numeric(frame[d+2])
		lab <- frame[d+3]
		if (nchar(lab) == 0) {
			highlwd <- 1
		} else {
			highlwd <- 3
		}		
		curvp <- viewport(x=unit(coords[1], "native"), y=unit(coords[2], "native"), width=margin[1], height=margin[2])
		grid.rect(x=unit(0.5, "snpc"), y=unit(0.5, "snpc"), width=unit(1.0, "snpc"), height=unit(1.0, "snpc"), gp=gpar(lwd=highlwd, col=rgb(0,0,0,alpha), fill=rgb(grey, grey, grey,alpha)), vp=curvp)
		grid.text(lab, y=unit(0.5, "snpc"), gp=gpar(fontsize=20, fontface="bold"), vp=curvp)
	}
	
	trianglepatch <- function(frame) {
		coords <- as.numeric(frame[1:d])
		maingrey <- as.numeric(frame[d+1])
		auxgrey <- as.numeric(frame[d+2])
		alpha <- as.numeric(frame[d+3])
		lab <- frame[d+4]
		if (nchar(lab) == 0) {
			highlwd <- 1
		} else {
			highlwd <- 3
		}
		curvp <- viewport(x=unit(coords[1], "native"), y=unit(coords[2], "native"), width=margin[1], height=margin[2])
		grid.polygon(x=unit(c(0,0,1), "snpc"), y=unit(c(1,0,0), "snpc"), gp=gpar(lwd=0, fill=rgb(maingrey, maingrey, maingrey, alpha)), vp=curvp)
		grid.polygon(x=unit(c(0,1,1), "snpc"), y=unit(c(1,1,0), "snpc"), gp=gpar(lwd=0, fill=rgb(auxgrey, auxgrey, auxgrey, alpha)), vp=curvp)	
		grid.rect(x=unit(0.5, "snpc"), y=unit(0.5, "snpc"), width=unit(1.0, "snpc"), height=unit(1.0, "snpc"), gp=gpar(lwd=highlwd, col=rgb(0,0,0,alpha), fill=rgb(1,1,1,0)), vp=curvp)	
		grid.text(lab, y=unit(0.5, "snpc"), gp=gpar(fontsize=20, fontface="bold"), vp=curvp)
	}
	
	# cbind handles NULL elements gracefully
	frame <- cbind(data, maingrey, auxgrey, alpha, labs)
	switch(method, 
		squarepatch=t(apply(frame, 1, squarepatch)),
		trianglepatch=t(apply(frame, 1, trianglepatch))
	)
		
	return(NULL)
}	
	
	
plotGmmOverlay <- function(gmm) {
	# plot gmm object, as defined in VBmix.
	# overlay it on existing grid graphic.
	# assumes the existence of an active grid graphic device,
	# with appropriately set scales
	# see mymvn2plot in VBmix package for reference
	if (!checkGmm(gmm)) {
		stop("argument should be a GMM object (see newGmm() in VBmix package for reference")
	}
	K <- length(gmm$w)
	nstep <- 20
	if(length(gmm$mean[[1]]) != 2) stop("only 2D GMMs can be used")
	
	# fill ellipses as a function of gmm$w
	for(i in 1:K) {
    	ev <- eigen(gmm$cov[[i]], symmetric = TRUE)
    	s <- sqrt(rev(sort(ev$values)))
    	V <- t(ev$vectors[, rev(order(ev$values))])
		theta <- (0:nstep) * (pi/(2 * nstep))
		x <- s[1] * cos(theta)
		y <- s[2] * sin(theta)    	
		xy <- cbind(c(x, -x, -x, x), c(y, y, -y, -y))
		xy <- xy %*% V
		xy <- sweep(xy, 2, gmm$mean[[i]], "+")
		# set everything in right order for drawing
		l <- length(x)
		L <- (l+1):(2*l)
		xy[L,] <- apply(xy[L,], 2, rev)
	
		L <- (3*l+1):(4*l)
		xy[L,] <- apply(xy[L,], 2, rev)
		w <- 1 - gmm$w[i]
		# log transform to weaken differences between low and high shades
		w <- log(1+w)
		grid.polygon(x=unit(xy[,1], "native"), y=unit(xy[,2], "native"), gp=gpar(lwd=2, fill=rgb(0, 0, 0, w)))
	}
	return(NULL)
}
	

checkGmm <- function(obj) {
	# check GMM structure
	# should be a list
	if (!is.list(obj)) return(FALSE)
	if (is.null(obj$w) || is.null(obj$mean) || is.null(obj$cov)) return(FALSE)
	if (!is.numeric(obj$w) || !is.list(obj$mean) || !is.list(obj$cov)) return(FALSE)
	return(TRUE)
}


getInterRGB <- function(vals, zeroColor, oneColor) {
	if (any(vals < 0 || vals > 1)) {
		stop("vals values should range in [0,1]")
	}

	labzero <- RGB2Lab(col2rgb(zeroColor) / 255)
	labone <- RGB2Lab(col2rgb(oneColor) / 255)
	
	interL <- vals * labone[1] + (1-vals) * labzero[1]
	intera <- vals * labone[2] + (1-vals) * labzero[2]
	interb <- vals * labone[3] + (1-vals) * labzero[3]
	
	rgbmat <- Lab2RGB(cbind(interL, intera, interb))
	rgbres <- rgb(rgbmat[,1], rgbmat[,2], rgbmat[,3])
	return(rgbres)
}



# plot a data.frame with color cells.
# data.frame values must be in [0,1] : the associated colors is then taken
# on the gradient from zeroColor to oneColor
# inspired from myImagePlot (http://www.phaget4.org/R/image_matrix.html)
colorPlot <- function(x, zeroColor, oneColor){
	if (any(x > 1) || any(x < 0)) {
		stop("x must take values in [0,1]")
	}
    yLabels <- rownames(x)
    xLabels <- colnames(x)
    title <-c()

	layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))

 	ColorRamp <- getInterRGB(seq(0,1,length=256), zeroColor, oneColor)
 	ColorLevels <- seq(0, 1, length=length(ColorRamp))

 	# Reverse Y axis
 	reverse <- nrow(x) : 1
 	yLabels <- yLabels[reverse]
 	x <- x[reverse,]

 	# Data Map
 	par(mar = c(3,5,2.5,2))
 	image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
 		ylab="", axes=FALSE, zlim=c(0, 1))

	axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=1.2)
 	axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
 		cex.axis=1.2)

 	# Color Scale
 	par(mar = c(3,2.5,2.5,2))
 	image(1, ColorLevels,
      	matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      	col=ColorRamp,
      	xlab="",ylab="",
      	xaxt="n")

 	layout(1)
 	return(NULL)
}




