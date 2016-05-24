setGeneric("corridor", function(x,speedProp=.75, circProp=.25, plot=FALSE, ...){standardGeneric("corridor")})
setMethod(f = "corridor",
	  signature=c(x=".MoveTrackSingle"),
	  definition=function(x, speedProp, circProp, plot, ...){
		  if (!isLonLat(x)) stop("Convert your dataset to longlat projection (use spTransform).")
		  if (n.locs(x)>2){
			  segLength <- apply(cbind(coordinates(x)[-n.locs(x),],coordinates(x)[-1,]), 1, function(y) spDistsN1(as.matrix(t(y[1:2])), as.matrix(t(y[3:4])), longlat=T)) ##kilometer
		  } else {stop("The data-set has less than 2 fixes")}

		  speed <- speed(x) #meter/second
		  speedQuant <- speed>=quantile(speed, probs=speedProp, na.rm=T)

		  segMid <- midPoint(coordinates(x)[-n.locs(x),],coordinates(x)[-1,])  
		  #Bart mid points works on spatialpoints* so we can omit the coordinates function here
		  ##Marco that is right, but how can you hand over the coordinates once cut at the end and once at the beginning
		  segRadius <- segLength/2
pkgLoad<-requireNamespace("maptools", quietly = TRUE) & requireNamespace("circular", quietly = TRUE)
		  if (pkgLoad ) {
			  tAzimuth <- maptools::trackAzimuth(coordinates(x))
			  pAzimuth <- ((180+tAzimuth)*2)%%360

			  inCircle <- lapply(1:(n.locs(x)-1),  function(i,segRadius, segMid){ which(spDistsN1(pts=as.matrix(segMid), pt=segMid[i,], longlat=T)<=segRadius[i])}, segMid=segMid, segRadius=segRadius) 

			  inpAzimuth <- lapply(inCircle, function(i, pAzimuth) pAzimuth[i], pAzimuth=pAzimuth)
			  circVar <- lapply(lapply(inpAzimuth, circular::circular, units="degrees"), circular::var.circular)
		  } else {
		    stop("Both the packages maptools and circular need to be installed")
		    }

		  if (length(unique(circVar))==1) {
		    
			  stop('There were less than the required 2 midpoints within the buffer along the whole track to calculate a variance.')
		  }else{circVar <- unlist(circVar)}

		  circVarQuant <- (circVar<=as.numeric(quantile(circVar[circVar!=0], probs=circProp, na.rm=T))) & circVar!=0
		  corrBehav <- which(speedQuant&circVarQuant)
		  inCircleCorrBehav  <- lapply(inCircle, function(v,behav) v[v%in%behav], behav=corrBehav)

		  nBehav <- unlist(lapply(inCircleCorrBehav, length))
		  nPoint <- unlist(lapply(inCircle, length))
		  corrPointsTmp <- which(nPoint>2 & (nBehav>(nPoint-nBehav)))

		  if(length(corrPointsTmp)==0) warning("No corridor points found!")
		  if(plot){

		    plot(coordinates(x), col="darkgrey", type="l")
		    if(length(corrPointsTmp)!=0){
		      corrPoints <- segMid[corrPointsTmp,, drop=F]
		      maxPoints <- max(nBehav[corrPointsTmp])
		      corrPointsCol <- apply(data.frame(corrPointsTmp), 1,function(i, tmp){grDevices::rgb(length(unlist(tmp[i]))/maxPoints, 1-length(unlist(tmp[i]))/maxPoints,1)}, tmp=inCircleCorrBehav)
		      
		      points(corrPoints[,1], corrPoints[,2], col=corrPointsCol, pch=20, ...)
		    }}
		  ##create a MoveBurst object, thus all information (speed, azimuth, ...) are stored in a way, that the actual segment midpoint is represented by the first coordinate of the segement
		  x$segMid_x <- c(segMid[,1], NA) #not elegant to add a zero again
		  x$segMid_y <- c(segMid[,2], NA) 
		  x$speed <- c(speed, NA)
		  x$azimuth <- c(tAzimuth, NA)
		  x$pseudoAzimuth <- c(pAzimuth, NA)
		  x$circVar <- c(circVar, NA)
		  corr <- rep("no corridor", n.locs(x))
		  corr[corrPointsTmp] <- "corridor"
		  return(move::burst(x, f=factor(corr[-length(corr)])))
	  })

setMethod(f = "corridor",
	  signature=c(x=".MoveTrackStack"),
	  definition=function(x, speedProp, circProp, ...){
		  lapply(split(x), corridor)
	  })
