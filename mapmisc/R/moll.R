

objFun = function(param, start, target){
	
	crsRot = CRS(paste(
					"+proj=ob_tran +o_proj=robin +o_lon_p=",
					param['lon'], " +o_lat_p=",
					param['lat'], 
					" +ellps=sphere +no_defs", sep=""))
	
	sum(
			(as.vector(spTransform(
							start, 
							crsRot				
					)@coords)*360/(2*pi)
			 - target
			)^2, na.rm=TRUE
	)
}


objFunAngle = function(param, start, target) {

	# start is a SpatialPoints object of length 2
	# target is x, y
	# point 1 should be due north of point 2 after transform
	
	crsRot = CRS(paste("+proj=ob_tran +o_proj=moll +o_lon_p=",
					param['lon'], " +o_lat_p=", param['lat'],
					" +lon_0=", param['wrap'], 
					" +lon_wrap=", param['wrap'],
					" +ellps=WGS84 +datum=WGS84 ",
					"+units=m +no_defs +towgs84=0,0,0",
					sep='')
	)
	
	startT = spTransform(
					start, 
				crsRot				
	)
	
	resDist = sum(
			(as.vector(startT@coords[1,] - target[c('x','y')])/100000
						)^2
	)
	
	resAngle = startT@coords[2,] - startT@coords[1,]
	resAngle = 90 - Arg(resAngle[1] + 1i * resAngle[2])*360/(2*pi) 
	
	resDist + resAngle^2
	
}

moll = function(x=0, angle=NULL, flip=FALSE) {
	
	
	if(is.numeric(x)){
		midX = x[1]
		if(length(x)==1) {
			midY = 0
		} else {
			midY = x[2]
		}
	} else {
		
		if(class(x)=='Extent'){
			xExtent = raster(x,crs=mapmisc::crsLL)
		} else if(length(grep("^SpatialPoints", class(x)))){
			if(length(x)==1){
				x = raster(extend(extent(x), 10^(-4)), crs=crs(x))
			}
		}
		xExtent = projectExtent(x, crsLL)
		midX = mean(c(xmin(xExtent),xmax(xExtent)))
		midY = mean(c(ymin(xExtent),ymax(xExtent)))
	}

	if(is.null(angle)){
		result = CRS(paste(
						"+proj=moll +lon_wrap=",
						midX, " +lon_0=",
						midX,
						" +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 ",
						"+units=m +no_defs +towgs84=0,0,0",
						sep=''))
	} else {

	newPole = geosphere::destPoint(
			c(midX, midY), b=angle, 
			d=pi/2, a=1, f=0
	)
	
#	stuff=geosphere::greatCircle(c(midX, midY), newPole, n=100)
	
#	geosphere::onGreatCircle(c(midX, midY), c(0,90), newPole)
	
	param = optim(
				c(lon=0,lat=0),
				objFun,
				start = SpatialPoints(
						newPole,
						proj4string=crsLL
						),
				target=c(NA,10^9)
				)
				
	param = param$par			
			
#	spTransform(
#			SpatialPoints(newPole, proj4string=crsLL),
#			CRS(paste(
#							"+proj=ob_tran +o_proj=longlat +o_lon_p=",
#							param['lon'], " +o_lat_p=",
#							param['lat'], 
#							" +lon_0=0 +ellps=WGS84 +no_defs", sep=""))
#	)@coords*360/(2*pi)
	
	
	newOrigin = spTransform(
			SpatialPoints(cbind(midX, midY), proj4string=crsLL),
			CRS(paste(
							"+proj=ob_tran +o_proj=longlat +o_lon_p=",
							param['lon'], " +o_lat_p=",
							param['lat'], 
							" +lon_0=0 +ellps=WGS84 +no_defs", sep=""))
	)@coords*360/(2*pi)

	# optimize again, with angle

	paramAgain = c(param, wrap=as.numeric(newOrigin[1]))

	start = cbind(midX, midY)
	
	start = SpatialPoints(rbind(
					start, 
					geosphere::destPoint(start, angle, d=1000)
					), proj4string=crsLL)
			
	newParam = optim(paramAgain, objFunAngle, start=start, 
			target=c(x=0, y=0, angle=as.numeric(angle))
	)$par		
		
	result = CRS(paste("+proj=ob_tran +o_proj=moll +o_lon_p=",
						newParam['lon'], " +o_lat_p=", newParam['lat'],
						" +lon_0=", newParam['wrap'], 
						" +lon_wrap=", newParam['wrap'],
						" +ellps=WGS84 +datum=WGS84 ",
						"+units=m +no_defs +towgs84=0,0,0",
						sep=''))
	
	}
	
	theBox = llCropBox(
			crs=result, res=1)

	if(is.character(flip)) {
		result = CRS(paste(as.character(result), " +axis=", flip, sep=''))
	} else {
		if(flip){
			if(midX < 0) {
				result = CRS(paste(as.character(result), "+axis=seu"))
			} else {
				result = CRS(paste(as.character(result), "+axis=wsu"))
			}
		} 
	}
	attributes(result)$crop = theBox$crop
	attributes(result)$ellipse = theBox$poly
	
	result
}

ocea = function(x, angle=0, flip=FALSE) {
	
	if(!requireNamespace('geosphere', quietly=TRUE) | 
			!requireNamespace('rgdal', quietly=TRUE)) {	
		warning("install geosphere and rgdal and to use ocea")
	}	
	
	
	northShift=0; eastShift=0; twistShift=0
	northShiftS = -60*60*northShift
	eastShiftS = 60*60*eastShift
	twistShiftS = 60*60*twistShift
	
	if(any(c(northShiftS, eastShiftS, twistShiftS)!= 0)){
		crsSphere	= CRS(paste(
						"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0,",
						northShiftS, ",",
						twistShiftS, ",",
						eastShiftS, ",0",
						sep=''))
	} else {
		crsSphere = mapmisc::crsLL
	}
	
	if(is.numeric(x)){
		if(length(x)==2) x = c(x,x+10^(-4))
		x = extent(x)
	}
	if(class(x)=='Extent'){
		xExtent = projectExtent(raster(x, crs=mapmisc::crsLL), crsSphere)
	} else {
		if(length(grep("^SpatialPoints", class(x)))){
			if(length(x)==1){
				x = raster(extend(extent(x), 10^(-4)), crs=crs(x))
			}
		}
		xExtent = projectExtent(x, crsSphere)
	}
	
	midY = mean(c(ymin(xExtent),ymax(xExtent)))
	midX = mean(c(xmin(xExtent),xmax(xExtent)))
	
	myCircle = geosphere::greatCircleBearing(
			c(midX, midY),
			angle, n=5
	)
	
	myEquator = geosphere::gcLon(myCircle[2,], myCircle[3,], 0)
	
	angleIntersection=geosphere::finalBearing(
			myCircle[2,], cbind(as.vector(myEquator),0)
	)
	
	whichPos = which.max(angleIntersection)
	
	myCrs = CRS(paste("+proj=ocea",
					" +lonc=", myEquator[whichPos], 
					" +alpha=", 180-angleIntersection[whichPos], 
					" +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84", 
					" +units=m +no_defs",
					" +towgs84=0,0,0,",  
					northShiftS, ",",
					twistShiftS, ",",
					eastShiftS, ",0",
					sep='') 
	)
	
	if(flip){
		if(midY < 0) {
			myCrs = CRS(paste(as.character(myCrs), "+axis=seu"))
		} else {
			myCrs = CRS(paste(as.character(myCrs), "+axis=nwu"))
		}
	}
	
	attributes(myCrs)$crop = llCropBox(
			crs=myCrs,
			keepInner=FALSE)$crop
	
	circleLLp = SpatialPoints(
			geosphere::greatCircle(
					myCircle[2,],
					myCircle[3,], n=500, sp=FALSE
			), proj4string=crsSphere)
	
	attributes(myCrs)$circleLL = spTransform(
			circleLLp, mapmisc::crsLL
	)
	
	attributes(myCrs)$circleTrans = spTransform(
			circleLLp, myCrs)
	
	myCrs
}
