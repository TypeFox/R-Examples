setGeneric("windowApply", function(x, FUN, SUMFUN, windowSize, ...) {
	   standardGeneric("windowApply")
})
setMethod("windowApply", signature(x = ".MoveTrackSingle", FUN = "function", 
				   SUMFUN = "function", windowSize = "numeric"), 
	  function(x, FUN, SUMFUN, windowSize, 
		   segmentWise = T, cluster = NULL, ...) {
		  if(n.locs(x)<windowSize)
			  stop("The window size can't be larger then the number of locations in move object")
		  # windowSize is the number of segments a window is long+1
		  RUNFUN <- function(x, windowSize, obj, WINFUN, segmentWise, ...) {
			  #		  message(x)# progress statement for debuging
			  res <- cbind(s = x:(x + windowSize - 1), WINFUN(obj[x:(x + windowSize - 1), 
									  ], ...))
	if (segmentWise) 
		res[nrow(res), ] <- NA
	return(res)
		  }
		  windowStarts <- 1:(n.locs(x) - windowSize + 1)
		  if (any(class(cluster) == "cluster")) {
			  res <- parLapply(cluster, as.list(windowStarts), FUN = RUNFUN, windowSize = windowSize, 
					   obj = x, WINFUN = FUN, segmentWise, ...)
		  } else {
			  res <- lapply(as.list(windowStarts), FUN = RUNFUN, windowSize = windowSize, 
					obj = x, WINFUN = FUN, segmentWise, ...)
		  }
		  res <- do.call("rbind", res)
		  res <- res[!is.na(res[, 2]), ]
		  val <- (aggregate(res[, -1], list(segNumber = res[, "s"]), SUMFUN))
		  val$nEstim <- tapply(res[,"s"], res[,"s"], length)[as.character(val$segNumber)]
		  return(val)
	  })

setGeneric("deltaParaOrth", function(mu, directionPoint, point) {
	   standardGeneric("deltaParaOrth")
	  })
setMethod("deltaParaOrth", signature(mu = "ANY", directionPoint = "numeric", 
				     point = "matrix"), function(mu, directionPoint, point) {
	directionPoint <- matrix(directionPoint, ncol = ncol(point), nrow = nrow(point), 
				 byrow = T)
	callGeneric()
	  })
setMethod("deltaParaOrth", signature(mu = "numeric", directionPoint = "matrix", 
				     point = "matrix"), function(mu, directionPoint, point) {
	mu <- matrix(mu, ncol = ncol(point), nrow = nrow(point), byrow = T)
	callGeneric()
	  })
setMethod("deltaParaOrth", signature(mu = "matrix", directionPoint = "matrix", 
				     point = "matrix"), function(mu, directionPoint, point) {
	C <- sqrt(rowSums((mu - directionPoint)^2))
	B <- sqrt(rowSums((point - directionPoint)^2))
	A <- sqrt(rowSums((point - mu)^2))
	tmp = ((A^2 + C^2-B^2)/(2 * A * C))
	sameLoc <- (rowSums(mu == point) == ncol(mu))
	sameLocDirection <- (rowSums(mu == directionPoint) == ncol(mu))
	tmp[sameLocDirection]<-sqrt(.5)
	if(any(sameLocDirection))
		warning("Brownian motion assumed, because no direction could be calculated")
	tmp[sameLoc] <- 0
	if (any(tmp > 1)) 
		tmp[tmp > 1] <- 1  # to prevent acosine from failing when above 1 due to float errors
	if (any(tmp < (-1))) 
		tmp[tmp < (-1)] <- (-1)  # to prevent acosine from failing when above 1 due to float errors
	theta1 <- acos(tmp)
	as.matrix(data.frame(deltaPara = A * cos(theta1),deltaOrth = A * sin(theta1)))

	  })

setGeneric("BGBvar", function(move, sdPara, sdOrth, locErr) {
	   standardGeneric("BGBvar")
	  })
setMethod("BGBvar", signature(move = ".MoveTrackSingle", sdPara = "numeric", 
			      sdOrth = "numeric", locErr = "numeric"), function(move, sdPara, sdOrth, locErr) {
	if ((n.locs(move)%%2) != 1) 
		stop("not an odd number of locations to BGBvar")
	g <- expand.grid(sdPara = sdPara, sdOrth = sdOrth)
	g$res <- NA
	t <- as.numeric(move@timestamps)/60#timestamps function
	if(length(locErr)==1)
		locErr <- rep(locErr, n.locs(move))
	is <- (1:n.locs(move))[(1:n.locs(move))%%2 == 0]
	alphas <- (t[is] - t[is - 1])/(t[is + 1] - t[is - 1])
	mus <- coordinates(move)[is - 1, ] + alphas * (coordinates(move)[is + 1, ] - 
						       coordinates(move)[is - 1, ])
	paraOrth <- deltaParaOrth(mus, coordinates(move)[is + 1, ], coordinates(move)[is, 
				  ])
    for (ii in 1:nrow(g)) {
	    sdPara <- g[ii, ]$sdPara
	    sdOrth <- g[ii, ]$sdOrth
	    sdParaTmp <- sqrt(alphas^2 * locErr[is + 1]^2 + (1 - alphas)^2 * locErr[is - 1]^2 + sdPara^2 * 
			      alphas * (1 - alphas) * (t[is + 1] - t[is - 1]))  # check if sdPara needs to be squared
	    sdOrthTmp <- sqrt(alphas^2 * locErr[is + 1]^2 + (1 - alphas)^2 * locErr[is - 1]^2 + sdOrth^2 * 
			      alphas * (1 - alphas) * (t[is + 1] - t[is - 1]))
	    g[ii, ]$res <- llBGBvar(cbind(sdParaTmp, sdOrthTmp)^2, paraOrth)
    }
    return(g)
	  })
setMethod("BGBvar", signature(move = ".MoveTrackSingle", sdPara = "missing", 
			      sdOrth = "missing", locErr = "numeric"), function(move, sdPara, sdOrth, locErr) {
	tmp <- optim(c(1, 1), function(x, move, locErr) {
		     -1 * BGBvar(move, x[1], x[2], locErr)$res
			      }, move = move, locErr = locErr, method = "L-BFGS-B", lower = 0, upper = 10^10)
	return(BGBvar(move, locErr = locErr, tmp$par[1], tmp$par[2]))
	  })
setGeneric("llBGBvar", function(sigm, paraOrth) {
	   # sigm is a matrix of variance values
	   standardGeneric("llBGBvar")
	  })
setMethod("llBGBvar", 
	  signature(sigm = "matrix", paraOrth = "matrix"), 
	  function(sigm, paraOrth) {
		  .Call('llBGBvar',sigm, paraOrth)
	  })
# setMethod('llBGBvar', signature(sigm = 'matrix', paraOrth = 'matrix'),
# function(sigm, paraOrth) { return(sum(unlist(lapply(1:nrow(sigm), function(i,
# sigm, paraOrth) { -log(2 * pi) - 0.5 * log(det(diag(2) * sigm[i, ])) - 0.5 *
# paraOrth[i, ] %*% solve(diag(2) * sigm[i, ]) %*% paraOrth[i, ] }, sigm =
# sigm, paraOrth = paraOrth)))) }) # this implementation is very slow about 20*
# slower than below there are float size differences between them but no idea
# what is better, differences are smaller than e-04

setGeneric("llBGBvarbreak", 
	   function(paraBreak, orthBreak,...) {
		   standardGeneric("llBGBvarbreak")
	   })
setMethod("llBGBvarbreak", 
	  signature(paraBreak = "missing", orthBreak = "missing"), 
	  function(paraBreak, orthBreak, paraBefore, 
		   orthBefore, paraOrth, errs, sdMul) {
		  sdOrthTmp <- errs + orthBefore^2 * sdMul
		  sdParaTmp <- errs + paraBefore^2 * sdMul
		  return(llBGBvar(cbind(sdParaTmp, sdOrthTmp), paraOrth))
	  })
setMethod("llBGBvarbreak", 
	  signature(
		    paraBreak = "numeric", orthBreak = "missing"), 
	  function(paraBreak, orthBreak,paraBefore, orthBefore, paraAfter,  
		   paraOrth, errs, sdMul ) {
		  para <- c(rep(paraBefore, paraBreak), rep(paraAfter, length(errs) - paraBreak))
		  sdOrthTmp <- errs + orthBefore^2 * sdMul
		  sdParaTmp <- errs + para^2 * sdMul
		  return(llBGBvar(cbind(sdParaTmp, sdOrthTmp), paraOrth))
	  })
setMethod("llBGBvarbreak", 
	  signature(paraBreak = "numeric",orthBreak="numeric" ), 
	  function( paraBreak, orthBreak,paraBefore, orthBefore, paraAfter, orthAfter, 
		   paraOrth, errs, sdMul) {
		  para <- c(rep(paraBefore, paraBreak), rep(paraAfter, length(errs) - paraBreak))
		  orth <- c(rep(orthBefore, orthBreak), rep(orthAfter, length(errs) - orthBreak))
		  sdOrthTmp <- errs + orth^2 * sdMul
		  sdParaTmp <- errs + para^2 * sdMul
		  return(llBGBvar(cbind(sdParaTmp, sdOrthTmp), paraOrth))
	  })
setMethod("llBGBvarbreak", 
	  signature(
		    paraBreak = "missing", orthBreak = "numeric"), 
	  function( paraBreak, orthBreak,paraBefore, orthBefore,  orthAfter, paraOrth, errs, sdMul) {
		  orth <- c(rep(orthBefore, orthBreak), rep(orthAfter, length(errs) - orthBreak))
		  sdOrthTmp <- errs + orth^2 * sdMul
		  sdParaTmp <- errs + paraBefore^2 * sdMul
		  return(llBGBvar(cbind(sdParaTmp, sdOrthTmp), paraOrth))
	  })
setGeneric("BGBvarbreak", 
	   function(move, locErr, margin, paraBreaks, orthBreaks, ...) {
		   standardGeneric("BGBvarbreak")
	   })
setMethod("BGBvarbreak", 
	  signature(move = ".MoveTrackSingle", locErr = "numeric", 
		    margin = "missing", paraBreaks = "numeric", orthBreaks = "numeric"), 
	  function(move, locErr, margin, paraBreaks, orthBreaks, ...) {
		  breaks <- expand.grid(paraBreaks = paraBreaks, orthBreaks = orthBreaks)
		  t <- as.numeric(move@timestamps)/60
		  if(length(locErr)==1)
			  locErr <- rep(locErr, n.locs(move))
		  is <- (1:n.locs(move))[(1:n.locs(move))%%2 == 0]
		  alphas <- (t[is] - t[is - 1])/(t[is + 1] - t[is - 1])
		  mus <- coordinates(move)[is - 1, ] + alphas * (coordinates(move)[is + 1, ] - 
								 coordinates(move)[is - 1, ])
		  paraOrth <- deltaParaOrth(mus, 
					    coordinates(move)[is + 1, ], 
					    coordinates(move)[is, ])
		  breaks$res <- NA
		  errs <- alphas^2 * locErr[is + 1]^2 + (1 - alphas)^2 * locErr[is - 1]^2
		  sdMul <- alphas * (1 - alphas) * (t[is + 1] - t[is - 1])
		  optims <- list()
		  for (i in 1:nrow(breaks)) {
			  paraBreak <- breaks[i, "paraBreaks"]
			  orthBreak <- breaks[i, "orthBreaks"]
			  init <- c(paraBefore = 100, orthBefore = 100)
			  otherArgs <- list(errs = errs, sdMul = sdMul, paraOrth = paraOrth)
			  s<-c(paraBreak='missing', orthBreak='missing')
			  if (!is.na(paraBreak)) {
				  init <- c(init, paraAfter = 100)
				  otherArgs <- c(otherArgs, list(paraBreak = floor(paraBreak/2)))
				  s['paraBreak']<-'numeric'
			  }
			  if (!is.na(orthBreak)) {
				  init <- c(init, orthAfter = 100)
				  otherArgs <- c(otherArgs, list(orthBreak = floor(orthBreak/2)))
				  s['orthBreak']<-'numeric'
			  }
			  # see if maybe it is a lot quicker to first find a method and then optimize so not every call needs to go through the standard generic
			  breaks[i, "res"] <- (tmp <- optim(init, function(x, otherArgs, ...) {
							    do.call(..., c(x, otherArgs))
									  }, 
									  method = "L-BFGS-B", 
									  control = list(fnscale = -1), 
									  lower = 0, 
									  upper = 10^10, 
									  otherArgs = otherArgs, what=getMethod('llBGBvarbreak', s)))$value
			  optims[[i]] <- tmp
		  }
		  breaks$BIC <- -2 * breaks$res + 
		  (2 + rowSums(!is.na(breaks[, c("paraBreaks", "orthBreaks")]))) * log(n.locs(move))  
		  brk <- breaks[which.min(breaks$BIC), ]
		  opt <- optims[[which.min(breaks$BIC)]]
		  if (any(opt$par == 0)) 
			  warning("Optimized to zero")
		  res <- cbind(paraSd = rep(opt$par["paraBefore"], n.locs(move)), orthSd = rep(opt$par["orthBefore"], 
											       n.locs(move)))
		  if (!is.na(brk$paraBreak)) {
			  res[1:nrow(res) >= brk$paraBreak, "paraSd"] <- opt$par["paraAfter"]
		  }
		  if (!is.na(brk$orthBreak)) {
			  res[1:nrow(res) >= brk$orthBreak, "orthSd"] <- opt$par["orthAfter"]
		  }
		  res[1:nrow(res) < min(c(paraBreaks, orthBreaks), na.rm = T), ] <- NA
		  res[1:nrow(res) >= max(c(paraBreaks, orthBreaks), na.rm = T), ] <- NA
		  rownames(res) <- NULL
		  return(as.data.frame(res))
	  })


setMethod("BGBvarbreak", 
	  signature(move = ".MoveTrackSingle", locErr = "numeric", 
		    margin = "numeric", paraBreaks = "missing", orthBreaks = "missing"), function(move, 
		    locErr, margin, paraBreaks, orthBreaks, ...) {
		  potentialBreaks <- 2:(n.locs(move) - 1)
		  marginUnevenPotentialBreaks <- potentialBreaks[potentialBreaks >= margin & potentialBreaks <= 
								 (1 + n.locs(move) - margin) & (potentialBreaks%%2) == 1]
		  orthBreaks <- c(NA, marginUnevenPotentialBreaks)
		  paraBreaks <- c(NA, marginUnevenPotentialBreaks)
		  callGeneric(move = move, locErr = locErr, paraBreaks = paraBreaks, orthBreaks = orthBreaks, 
			      ...)
	  })
setGeneric("BGB", function(move, raster, locErr) {
	   standardGeneric("BGB")
	  })
setMethod("BGB", 
	  signature(move = ".MoveTrackSingle", raster = "numeric", locErr = "numeric"), 
	  function(move, raster, locErr) {
		  xRange <- range(coordinates(move)[, 1])
		  yRange <- range(coordinates(move)[, 2])
		  # make the range a bit larger
		  xRange <- xRange + c(-0.2, 0.2) * diff(xRange)
		  yRange <- yRange + c(-0.2, 0.2) * diff(yRange)
		  # make the range fit the raster cells
		  xRange <- xRange + c(-0.5, 0.5) * (raster - diff(xRange)%%raster)
		  yRange <- yRange + c(-0.5, 0.5) * (raster - diff(yRange)%%raster)
		  ex <- extent(c(xRange, yRange))
		  raster <- raster(ncols = diff(xRange)%/%raster, nrows = diff(yRange)%/%raster, 
				   crs = proj4string(move), ex)
		  callGeneric()
	  })
setMethod("BGB", 
	  signature(move = ".MoveTrackSingle", raster = "RasterLayer", locErr = "numeric"), 
	  function(move, raster, locErr) {
		  time.step <- min(timeLag(move, units="mins"))/15.123
		  x <- coordinates(move)[, 1]
		  y <- coordinates(move)[, 2]
		  location.error <- rep(locErr, length(x))
		  t <- as.numeric(move@timestamps)/60
		  tm <- t[1] + time.step/2
		  g <- BGBvar(move, locErr = locErr)
		  sigmaSeg_orth <- rep(g$sdOrth, length(x))
		  sigmaSeg_para <- rep(g$sdPara, length(x))
		  int <- 0
		  for (i in 1:(length(x) - 1)) {
			  # tm <- 0 #
			  while (tm <= t[i + 1]) {
				  alpha <- (tm - t[i])/(t[i + 1] - t[i])
				  mu.x <- x[i] + alpha * (x[i + 1] - x[i])
				  mu.y <- y[i] + alpha * (y[i + 1] - y[i])
				  paraOrth <- deltaParaOrth(c(mu.x, mu.y), c(x[i + 1], y[i + 1]), coordinates(raster))
				  sigma_para <- sqrt((t[i + 1] - t[i]) * alpha * (1 - alpha) * (sigmaSeg_para[i])^2 + 
						     ((1 - alpha)^2) * (location.error[i]^2) + (alpha^2) * (location.error[i + 
													    1]^2))
		sigma_orth <- sqrt((t[i + 1] - t[i]) * alpha * (1 - alpha) * (sigmaSeg_orth[i])^2 + 
				   ((1 - alpha)^2) * (location.error[i]^2) + (alpha^2) * (location.error[i + 
											  1]^2))
		if (any(is.nan(c(sigma_orth, sigma_para)))) 
			browser()
		out = (1/(2 * pi * sigma_para * sigma_orth)) * exp(-0.5 * ((paraOrth$para^2/sigma_para^2) + 
									   (paraOrth$orth^2/sigma_orth^2)))
		if (any(is.nan(out))) 
			stop('NAN in output')
		int <- int + out/sum(out)
		if (any(is.nan(int))) 
			stop('NAN in output')
		tm <- tm + time.step
			  }
		  }
		  # Scaling probabilities so they sum to 1.0
		  values(raster) <- int/sum(int)
		  return(raster)
	  })

setMethod("[", "dBGBvariance", function(x, i, j, ..., drop = TRUE) {
	  x@paraSd <- x@paraSd[i]
	  x@orthSd <- x@orthSd[i]
	  x@nEstim <- x@nEstim[i]
	  x@segInterest <- x@segInterest[i]
	  callNextMethod()
	  })


setGeneric("dynBGBvariance", 
	   function(move, locErr, margin, windowSize, ...) {
		   standardGeneric("dynBGBvariance")
	   })
setMethod("dynBGBvariance", signature(move = ".MoveTrackSingle", locErr = "numeric", 
				      margin = "numeric", windowSize = "numeric"), function(move, locErr, margin, windowSize, 
				      ...) {
	vars <- windowApply(move, BGBvarbreak, function(x){sqrt(mean(x^2))}, locErr = locErr, windowSize = windowSize, # sum to average variances
			    ..., margin = margin) 
	tmp <- as.data.frame(matrix(NA, ncol = 3, nrow = n.locs(move)))

	tmp[vars$segNumber, ] <- (vars[, 2:4])
	names(tmp) <- names(vars[, 2:4])
	new("dBGBvariance", 
	    move, 
	    margin = margin,
	    windowSize = windowSize,
	    paraSd = tmp$para, 
	    orthSd = tmp$orth, 
	    nEstim = tmp$nEstim, 
	    segInterest = !(tmp$nEstim != max(tmp$nEstim, na.rm = T) | is.na(tmp$nEstim)))
	   })
setGeneric("dynBGB", 
	   function(move, raster, locErr, ...) {
		   standardGeneric("dynBGB")
	   })
setMethod("dynBGB", signature(move = ".MoveTrackSingle", raster = "RasterLayer", 
			      locErr = "numeric"), 
	  function(move, raster, locErr, margin, windowSize, ...) {
		  move <- dynBGBvariance(move = move, locErr = locErr, margin=margin, windowSize=windowSize,...)
		  callGeneric(move = move, raster = raster, locErr = locErr,...)
	  })

setMethod("dynBGB", signature(move = ".MoveTrackSingle", raster = "ANY", 
			      locErr = "character"), 
	  function(move, raster, locErr, ...) {
		  locErr <- do.call("$", list(move, locErr))
		  if(is.null(locErr))
			  stop('column indicated for locErr probably does not exist')
		  callGeneric(move = move, raster = raster, locErr = locErr,...)
	  })

setMethod("dynBGB", signature(move = ".MoveTrackSingle", raster = "numeric", locErr = "ANY" 
			      ), function(move, raster, 
			      locErr, ext,...) {
	r<-raster
	raster <- raster(extent(.extcalc(obj = move, ext = ext)))
	res(raster)<-r
	proj4string(raster)<- proj4string(move)
	callGeneric(move = move, raster = raster, locErr = locErr,...)
	  })

setMethod("dynBGB", signature(move = ".MoveTrackSingle", raster = "missing", locErr = "ANY" 
			      ), function(move, raster, 
			      locErr, dimSize,ext,...) {
	raster <- raster::raster(extent(.extcalc(obj = move, ext = ext)))
	res(raster)<-max(diff(t(bbox(raster))))/dimSize
	proj4string(raster)<- proj4string(move)
	callGeneric(move = move, raster = raster, locErr = locErr,...)
	  })

setMethod("dynBGB", 
	  signature(move = "dBGBvariance", raster = "RasterLayer", 
		    locErr = "numeric"), 
	  function(move, raster, locErr, timeStep, ...) {
		  if(isLonLat(move)) stop("You can not use longitude latitude projection for this function. To transform your coordinates use the spTransform function. \n")
		  if(!equalProj(list(raster,move))) #check equal projection of raster and Move
			  stop(paste("The projection of the raster and the Move object are not equal. \n raster:", proj4string(raster), "\n object:", proj4string(move), "\n"))
		  pointsInterest <- move@segInterest | rev(move@segInterest)
		  t <- as.numeric(move@timestamps)/60
		  if(missing(timeStep))
			  timeStep<-min(diff(t))/20.1
		  #dyn.load("bbmm.so")
		  values(raster) <- .Call("bgb", 
					  coordinates(move)[pointsInterest, 1], 
					  coordinates(move)[pointsInterest, 2], 
					  move@paraSd[pointsInterest], 
					  move@orthSd[pointsInterest], 
					  t[pointsInterest], 
					  rep(locErr, sum(pointsInterest)), 
					  unique(coordinates(raster)[, 1]), 
					  sort(unique(coordinates(raster)[, 2])), 
					  timeStep, 
					  5# number of sd intergration distance
					  )
		  raster<-raster/cellStats(raster, sum)
		  res<-new("dynBGB", raster, var=move, method="dynBGB")
		  return(res)
	  })
#setGeneric("BGBUD", function(move, raster, locErr, maxInt) {
#	   standardGeneric("BGBUD")
#	  })
#setMethod("BGBUD", signature(move = "dBGBvariance", raster = "numeric", 
#			     locErr = "numeric", maxInt = "numeric"), function(move, raster, locErr, maxInt) {
#	xRange <- range(coordinates(move)[, 1])
#	yRange <- range(coordinates(move)[, 2])
#	# make the range a bit larger
#	xRange <- xRange + c(-1.9, 1.9) * diff(xRange)
#	yRange <- yRange + c(-1.9, 1.9) * diff(yRange)
#	# make the range fit the raster cells
#	xRange <- xRange + c(-0.5, 0.5) * (raster - diff(xRange)%%raster)
#	yRange <- yRange + c(-0.5, 0.5) * (raster - diff(yRange)%%raster)
#	ex <- extent(c(xRange, yRange))
#	raster <- raster(ncols = diff(xRange)%/%raster, nrows = diff(yRange)%/%raster, 
#			 crs = proj4string(move), ex)
#	callGeneric()
#	  })
#setMethod("BGBUD", signature(move = "dBGBvariance", raster = "RasterLayer", 
#			     locErr = "numeric", maxInt = "numeric"), function(move, raster, locErr, maxInt) {
#	dyn.load("bbmm.so")
#	pointsInterest <- move@segInterest | rev(move@segInterest)
#	t <- as.numeric(move@timestamps)/60
#	values(raster) <- .Call("bgb", coordinates(move)[pointsInterest, 1], coordinates(move)[pointsInterest, 
#				2], move@paraSd[pointsInterest], move@orthSd[pointsInterest], t, rep(locErr, 
#												     sum(pointsInterest)), unique(coordinates(raster)[, 1]), sort(unique(coordinates(raster)[, 
#																					 2])), min(diff(t))/90.1, maxInt,4)
#    return(raster)
#	  }) 
setAs( "dBGBvariance","data.frame", function(from){ return(data.frame(as(as(from,".MoveTrack"),"data.frame"), paraSd=from@paraSd, orthSd=from@orthSd, margin=from@margin, windowSize=from@windowSize))})

