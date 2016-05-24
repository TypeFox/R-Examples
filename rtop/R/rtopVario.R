rtopVariogram.SpatialPolygonsDataFrame = function(object, ... ) {
  if (missing(object))  stop("rtopVariogram: Observations are missing")
  obs = object@data
  coordinates(obs) = coordinates(object)
  if ("Shape_Area" %in% names(object)) {
    obs$area = object$Shape_Area
  } else obs$area = unlist(lapply(object@polygons,FUN = function(poly) poly@area))
  rtopVariogram(obs, ...)
}





rtopVariogram.SpatialPointsDataFrame = function(object, formulaString, params=list(), cloud, abins, dbins, ...) {
  # If params is intamapParams, they will here be included in rtopParams
  if (!inherits(params, "rtopParams"))  params = getRtopParams(params, ...)
  # amul refers to the number of areal bins per order of magnitude
  # dmul refers to the number of distance bins per order of magnitude
  amul = params$amul
  dmul = params$dmul
  if (missing(cloud)) cloud = params$cloud
  
  observations = object
  if (missing(observations)) stop("rtopVariogram: Observations are missing")
  if (!("area") %in% names(observations)) {
    if ("AREA" %in% names(observations)) {
      observations$area = observations$AREA
    } else {
      stop("rtopVariogram: Observations do not include the area of the polygons")
    }
  }
  #  stop("rtopVariogram: Observations do not include area (polygons) or length (lines)")
  if (missing(formulaString)) {
    if ("obs" %in% names(observations)) { 
      formulaString = "obs ~ 1" 
    } else if ("value" %in% names(observations)) {
      formulaString = "value ~ 1" 
    } else if (length(names(observations@data)) == 1) {
      formulaString = paste(names(observations@data),"~ 1")      
    } else stop("formulaString is missing and cannot be found from data")
    warning(paste("formulaString missing, using",formulaString))      
  }
  if (!inherits(formulaString,"formula")) formulaString = as.formula(formulaString)
  
  clvar = variogram(formulaString, observations, cloud = TRUE, ...)
  .BigInt = attr(clvar, ".BigInt")
  clvar$ord = clvar$np
  clvar = as.data.frame(clvar)
  
  clvar$a1 = observations@data$area[clvar$left]
  clvar$a2 = observations@data$area[clvar$right]
  
  
  if (cloud) {
    clvar$acl1 = clvar$left
    clvar$acl2 = clvar$right
    clvar = clvar[,-which(names(clvar) %in% c("left","right"))]
    var3d = clvar
    class(var3d) = c("rtopVariogramCloud","data.frame")
    attr(var3d, ".BigInt") = .BigInt
    var3d$np = 1
  } else {
    abins = adfunc(NULL, observations, amul)
    dbins = dfunc(NULL, observations, dmul)
    observations$acl = findInterval(observations$area, abins)
    clvar$acl1 = observations$acl[clvar$left]
    clvar$acl2 = observations$acl[clvar$right]
    
    ich = which(clvar$acl1 > clvar$acl2)
    acl1c = clvar$acl1
    clvar$acl1[ich] = clvar$acl2[ich]
    clvar$acl2[ich] = acl1c[ich]
    
    clvar$dbin = findInterval(clvar$dist, dbins)
    clvar$np = 1
    varnp = aggregate(list(np = clvar$np),list(acl1 = clvar$acl1,acl2 = clvar$acl2,dbin = clvar$dbin),sum)
    var3d = data.frame(np=varnp$np,aggregate(list(dist=clvar$dist,gamma=clvar$gamma,a1 = clvar$a1, a2=clvar$a2),
                                             list(acl1 = clvar$acl1,acl2 = clvar$acl2,dbin = clvar$dbin),mean))
    class(var3d) = c("rtopVariogram","data.frame")
  }
  var3d
}



# Alternative binning:
#  x <- matrix(rnorm(30000), ncol=3)
#  breaks <- seq(-1, 1, length=5)
#  xints <- data.frame(
#  x1=cut(x[, 1], breaks=breaks),
#  x2=cut(x[, 2], breaks=breaks),
#  x3=cut(x[, 3], breaks=breaks))
#  table(complete.cases(xints))
#  xtabs(~ ., xints)




###############################
rtopVariogram.rtop = function(object, params = list(), ... ) {
  params = getRtopParams(object$params, newPar = params, ...)
  observations = object$observations
  formulaString = object$formulaString
  
  #calling rtopVariogram.SpatialPolygonsDataFrame
  var3d = rtopVariogram(observations, formulaString, params,...)
  if (inherits(var3d, "rtopVariogramCloud")) object$variogramCloud = var3d else object$variogram = var3d
  object
}


################################


rtopVariogram.STSDF = function(object, formulaString, params = list(), cloud, abins, 
                               dbins, data.table = FALSE, ...) {
  if (!requireNamespace("spacetime")) stop("spacetime not available")
  if (!inherits(params, "rtopParams")) 
  params = getRtopParams(params, ...)
  amul = params$amul
  dmul = params$dmul
  if (missing(cloud)) 
    cloud = params$cloud
  observations = object
  debug.level = params$debug.level
  if (missing(observations)) 
    stop("rtopVariogram: Observations are missing")
  if (!("area") %in% names(observations@sp)) {
    if ("AREA" %in% names(observations@sp)) {
      observations@sp$area = observations@sp$AREA
    } else {
      stop("rtopVariogram: Observations do not include the area of the polygons")
    }
  }
  if (missing(formulaString)) {
    if ("obs" %in% names(observations@data)) {
      formulaString = "obs ~ 1"
    }
    else if ("value" %in% names(observations@data)) {
      formulaString = "value ~ 1"
    }
    else if (length(names(observations@data)) == 1) {
      formulaString = paste(names(observations@data), "~ 1")
    }
    else stop("formulaString is missing and cannot be found from data")
    warning(paste("formulaString missing, using", formulaString))
  }
  if (!inherits(formulaString, "formula")) 
    formulaString = as.formula(formulaString)
  depvar = as.character(formulaString[[2]])
  
  diffsN1 = function(x, y) y-x
  
  nspace = dim(observations)[1]
  observations@sp$vindex = 1:nspace
  ntime = dim(observations)[2]
  observations@sp$vindex = vindex = 1:nspace
  obsdf = as.data.frame(observations)[,c("timeIndex","vindex",depvar)]
  vmat = matrix(0, nrow = nspace, ncol = nspace)
  indmat = vmat
  if (interactive() & debug.level) pb <- txtProgressBar(1, ntime, style = 3)
  if (data.table && !requireNamespace("data.table")) {
    warning("data.table not available, continuing without")
    data.table = FALSE
  }
  if(data.table){
    message("Converting STSDF class to data.table class")
    observationsDT <- data.table::data.table(obsdf,key=c("timeIndex"))
    for (ind in 1:ntime) {
      ppq <- observationsDT[list(ind)]	
      nspace1 <-  dim(ppq)[1]
      ff <-  matrix(unlist(lapply(ppq[[depvar]],  FUN = function(x) diffsN1(x, ppq[[depvar]]))), 
                    byrow = TRUE, ncol = nspace1)
      ff <-  (ff^2)/2
      findx <-  ppq[,vindex]
      vmat[findx,findx] <-  vmat[findx,findx] + ff
      if (any(is.na(ff))) stop("na-values in covariance matrix")
      indmat[findx,findx] <-  indmat[findx,findx] + 1
      if (interactive() & debug.level) setTxtProgressBar(pb, ind)
    }
  }else{
    for (ind in 1:ntime) {
      ppq = obsdf[obsdf$timeIndex == ind,]
      nspace1 = dim(ppq)[1]
    ff = matrix(unlist(lapply(ppq[,depvar],  FUN = function(x) diffsN1(x, ppq[,depvar]))), 
                byrow = TRUE, ncol = nspace1)
    ff = (ff^2)/2
    findx = ppq$vindex
    vmat[findx,findx] = vmat[findx,findx] + ff
    if (any(is.na(ff))) stop("na-values in covariance matrix")
    indmat[findx,findx] = indmat[findx,findx] + 1
    if (interactive() & debug.level) setTxtProgressBar(pb, ind)
  }
  }
  if (interactive() & debug.level) close(pb)
  vmat = vmat/indmat
  obssp = observations@sp 
  dmat = spDists(obssp, obssp) 
  
  vario = matrix(NA, ncol = 7, nrow = nspace*(nspace-1)/2)
  icount = 0
  for (istat in 1:(nspace -1)) {
    njs = nspace-istat
    vario[(icount+1):(icount+njs),] = matrix(c(dmat[istat, (istat+1):nspace],
                                               vmat[istat, (istat+1):nspace], rep(obssp$area[istat], njs),
                                               obssp$area[(istat+1):nspace], rep(istat, njs), (istat+1):nspace,
                                               indmat[istat, (istat+1):nspace]), ncol = 7)
    icount = icount + njs
  } 
  vario = data.frame(vario)
  names(vario) = c("dist", "gamma", "a1", "a2", "acl1", "acl2", "np")
  vario = vario[!is.na(vario$gamma),]
  
  if (cloud) {
    var3d = vario
    class(var3d) = c("rtopVariogramCloud","data.frame")
  } else {
    abins = adfunc(NULL, obssp, amul)
    dbins = dfunc(NULL, obssp, dmul)
    obssp$acl = findInterval(obssp$area, abins)
    vario$acl1 = obssp$acl[vario$acl1]
    vario$acl2 = obssp$acl[vario$acl2]
    
    ich = which(vario$acl1 > vario$acl2)
    acl1c = vario$acl1
    vario$acl1[ich] = vario$acl2[ich]
    vario$acl2[ich] = acl1c[ich]
    
    vario$dbin = findInterval(vario$dist, dbins)
    varnp = aggregate(list(np = vario$np), list(acl1 = vario$acl1, acl2 = vario$acl2, dbin = vario$dbin),sum)
    var3d = data.frame(np = varnp$np, aggregate(list(dist = vario$dist, 
                                                     gamma=vario$gamma, a1 = vario$a1, a2 = vario$a2),
                                                list(acl1 = vario$acl1, acl2 = vario$acl2, dbin = vario$dbin),mean))
    class(var3d) = c("rtopVariogram","data.frame")
  }
  var3d
}






