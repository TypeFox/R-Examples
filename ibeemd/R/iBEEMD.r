# Created by Mao-Gui Hu (humg@lreis.ac.cn), 2013.10.
# 2d EMD for irregular geographical polygons
require(sp)
require(fields)

# Local min/max extremes
LocalExtremes <- function(spXXDf, nbList, boundary, bPlateauAndBeam=TRUE)
{
# spXXDf: SpatialPointsDataFrame or SpatialPolygonsDataFrame
# nbList: neighbours of each objects in spXXDf [list]
# boundary: indicate whether an object in spXXDf is located at boundary [vector]
# bPlateauAndBeam: whether plateau/pan and beam/channel should be treated as extremes

	if(!bPlateauAndBeam) {
		locMax <- c()
		locMin <- c()
		for(i in 1:nrow(spXXDf)) {
			if(boundary[i]) next
			
			curValue <- as.vector(spXXDf@data$value[i])
			nbValue <- as.vector(spXXDf@data$value[nbList[[i]]])
			if(all(curValue > nbValue)) locMax <- c(locMax, i)
			if(all(curValue < nbValue)) locMin <- c(locMin, i)
		}
	}
	else {
		nDiffWithMin <- rep(NA, nrow(spXXDf))
		nDiffWithMax <- rep(NA, nrow(spXXDf))
		for(i in 1:nrow(spXXDf)) {
			if(boundary[i]) next
			curV <- as.vector(spXXDf@data$value[i])
			minV <- min(c(curV,as.vector(spXXDf@data$value[nbList[[i]]])))
			maxV <- max(c(curV,as.vector(spXXDf@data$value[nbList[[i]]])))
			
			nDiffWithMin[i] <- curV - minV
			nDiffWithMax[i] <- curV - maxV
		}
		locMin <- which(nDiffWithMin == 0)
		locMax <- which(nDiffWithMax == 0)
	}
	return(list(locMax=locMax, locMin=locMin))
}

PlotFit <- function(xy, fitrslt, nx=64, ny=64) {
	par(mfrow=c(2,2))
	quilt.plot(xy, fitrslt$upperEnv, main="Upper envelope", nx=nx, ny=ny)
	quilt.plot(xy, fitrslt$lowerEnv, main="Lower envelope", nx=nx, ny=ny)
	quilt.plot(xy, fitrslt$meanSurf, main="Mean", nx=nx, ny=ny)
	quilt.plot(xy, fitrslt$residual, main="Residual", nx=nx, ny=ny)
}

# Spline interpolation of upper & lower envelopes
SurfaceFit <- function(spPointsDf, locExtr, fmodel="thinplate", fig=FALSE, nx=64, ny=64)
{
	# scale coordinates to an appropriate range
	ScaleCoord <- function(xy, nmax=100) {		
		return(xy / max(abs(xy)) * nmax)
	}
	
# 	cat("Surface fit.\n")
	LocalFit <- function()
	{
		if(fmodel == "cubic" || fmodel == "multiquadric" 
				|| fmodel == "gaussian" || fmodel == "thinplate") {
			xy <- ScaleCoord(spPointsDf@coords)
			fitMax <- rbfcreate(xy[locExtr$locMax,], spPointsDf@data$value[locExtr$locMax], rbffun=fmodel)
			fitMin <- rbfcreate(xy[locExtr$locMin,], spPointsDf@data$value[locExtr$locMin], rbffun=fmodel)
			upperEnv <- as.vector(rbfinterp(fitMax, xy))
			lowerEnv <- as.vector(rbfinterp(fitMin, xy))
		}
		else {
			stop("invalid surface fitting method.")
		}
		
		meanSurf <- (upperEnv + lowerEnv)/2
		residual <- spPointsDf@data$value - meanSurf
		
		fitrslt <- list(#fitMax=fitMax, fitMin=fitMin, 
					upperEnv=upperEnv, lowerEnv=lowerEnv, 
					meanSurf=meanSurf, residual=residual)
		
		if(fig) {
			PlotFit(spPointsDf@coords, fitrslt, nx=nx, ny=ny)
		}
		return(fitrslt)
	}
	
	tryCatch({
		fitrslt <- LocalFit()
		return(fitrslt)
	},
	error=function(x) {
		print(x) 
		return(NULL)
	}
	)
}

ExtractIMFi2D <- function(spPointsDf, nbList, boundary, tolSift=0.05, minExtrema=5,
		fmodel="thinplate", bPlateauAndBeam=TRUE, fig=FALSE, nx=64, ny=64)
{
# spPointsDf: SpatialPointsDataFrame
# nbList: neighbours of each objects in spXXDf [list]
# boundary: indicate whether an object in spXXDf is located at boundary [vector]
# bPlateauAndBeam: whether plateau/pan and beam/channel should be treated as extremes

# 	cat("Extract IMF.\n")	
	SD <- Inf
	nSift <- 0	# max sift iteration
	
	orgValue <- spPointsDf@data$value
	residual <- orgValue
	imf <- c()
	while(1) {
		spPointsDf@data$value <- residual
		locExtr <- LocalExtremes(spPointsDf, nbList, boundary, bPlateauAndBeam)
		if(length(locExtr$locMax) < minExtrema || length(locExtr$locMin) < minExtrema) break
		
		fit <- SurfaceFit(spPointsDf, locExtr, fmodel=fmodel, fig=fig, nx=nx, ny=ny)
		if(is.null(fit) || any(is.nan(fit$residual))) break
		SD <- sum((residual-fit$residual)^2)/sum(residual^2)
		
		nSift <- nSift+1
		if((SD < tolSift && nSift > 1) || nSift > 5) break
		
		imf <- fit$residual
		residual <- imf
	}
	if(is.null(imf)) 
		return(NULL)
	else
		residual <- orgValue-imf
	
	if(fig) {
		xy <- spPointsDf@coords
		par(mfrow=c(1,2))
		quilt.plot(xy, imf, main="IMF")
		quilt.plot(xy, residual, main="Residual")
	}
	
	return(list(imf=imf, residual=residual))
}

MeanIMF <- function(imfs)
{
	if(length(imfs) == 1) return(imfs[[1]])
	
	nimf <- rep(0, length(imfs))
	for(i in 1:length(imfs)) {
		nimf[i] <- NCOL(imfs[[i]]$imf)
	}
	
	uniqLen <- sort(unique(nimf))
	uniqLenNum <- c()
	for(i in uniqLen) {
		uniqLenNum <- c(uniqLenNum, sum(nimf==i))
	}
	selImfs <- imfs[nimf == uniqLen[which.max(uniqLenNum)]]
	mnImf <- 0
	mnTrend <- 0
	for(i in 1:length(selImfs)) {
		mnImf <- mnImf+selImfs[[i]]$imf
		mnTrend <- mnTrend+selImfs[[i]]$trend		
	}
	mnImf <- mnImf/length(selImfs)
	mnTrend <- mnTrend/length(selImfs)
	
	return(list(imf=mnImf, trend=mnTrend))
}

# decompostion of a SpatialPolygonsDataFrame by iBEEMD
# a matrix will be returned: 
#  the first column is the original value;
#  the last column is the decomposed residual (global tredn);
#  columns from the second to the last but one is decmompsed IMFs.
iBEEMD <- function(
		spPolysDf, 				# interested spatial data; it is a SpatialPolygonsDataFrame object.
		valueField = names(spPolysDf)[1],  	# interested field value in above spatial data.
		nMaxIMF = 10, 			# maximum number of IMFs (default is 10)
		tolSift = 0.05, 		# sift tolerence (default is 0.05)
		neemd = 1000,			# EEMD iteration number (default is 1000)
		wnsd = 0.05, 			# standard deviation of added noise; it is a ratio to the standard deviation of above data.
		fmodel = "thinplate", 	# surface fitting function ("thinplate", "gaussian", "cubic", "multiquadric")
		fig = TRUE				# show decomposed result as an image
		)
{
	stopifnot(is(spPolysDf, "SpatialPolygonsDataFrame"))
	spPolysDfBak <- spPolysDf
	
	minExtrema <- 5			# minimum number of nodes to fit a surface
	bPlateauAndBeam <- TRUE # whether plateau/pan and beam/channel should be treated as extremes
	spClipBnd <- NULL		# if the spPolysDf is a SpatialPointsDataFrame instead of SpatialPolygonsDataFrame,
							#  a SpatialPolygonsDataFrame dataset will be created by voronoi partition. Then, 
							#  spClipBnd is a SpatialPolygonsDataFrame containing only one polygon to define the boundary.
	
	# reserved for points
	bInputPt <- is(spPolysDf, "SpatialPointsDataFrame")
	if(bInputPt) {
		cat("Create polygons from points ...\n")
		spInputPt <- spPolysDf
		spPolysDf <- CreateVoronoi(spPolysDf, bDelaunay=FALSE, bClipVor=TRUE, spClipBnd=spClipBnd)
	}
	
	spPolysDf@data <- data.frame(id=1:length(spPolysDf), value=spPolysDf@data[,valueField])
	noiseSD <- sd(spPolysDf@data$value) * wnsd
	
	if(!bInputPt) {
		spPointsDf <- SPolygonsDfToSPointsDf(spPolysDf)
	} else {
		spPointsDf <- spInputPt
		spPointsDf@data <- spPolysDf@data
	}
	
	orgLen <- length(spPolysDf)
	orgValue <- spPolysDf@data$value
	spPolysDf <- createSPComment(spPolysDf)
	
	nbList <- NeighborPolys(spPolysDf)
	boundary <- LocatedAtBoundary(spPolysDf, bRemoveHole=TRUE)
	
	cat("EEMD iteration ...\n")
	
	if(noiseSD == 0 || neemd <= 0) neemd <- 1
	imfs <- list()
	bakValue <- spPointsDf@data$value
		
	pb <- txtProgressBar(min=0, max=neemd, width=50, style=3)
	for(k in 1:neemd)
	{
		# extract imfs
		imf <- c()
		trend <- c()
		spPointsDf@data$value <- bakValue+rnorm(length(spPointsDf),sd=noiseSD)
		for(i in 1:nMaxIMF)
		{
			extimf <- ExtractIMFi2D(spPointsDf, nbList, boundary, tolSift=tolSift, minExtrema=minExtrema,
 				bPlateauAndBeam=bPlateauAndBeam, fig=FALSE, fmodel=fmodel)
			if(is.null(extimf)) break
			
			imf <- cbind(imf, as.vector(extimf$imf))
			trend <- extimf$residual
			spPointsDf@data$value <- trend
			if(sum(trend) == 0) break
		}
		if(is.null(imf)) next
		imfs[[length(imfs)+1]] <- list(imf=imf, trend=as.vector(trend))
		setTxtProgressBar(pb, k)
	}
	close(pb)
	
	if(length(imfs) == 0) {
		warning("no imf decomposed.")
		return(NULL)
	}
	
	if(!is.null(imfs))
		mnImfTrend <- MeanIMF(imfs)
	
	imfDf <- as.data.frame(mnImfTrend$imf)
	colnames(imfDf) <- paste("imf", 1:ncol(imfDf), sep="")
	spPolysDfBak@data <- data.frame(original=orgValue, imfDf, trend=mnImfTrend$trend)
	names(spPolysDfBak@data)[1] <- valueField
	
	if(fig) {
		cat("Plotting ...\n")
		print(spplot(spPolysDfBak, col.regions = terrain.colors(256)))
	}
	
	return(spPolysDfBak)
}
