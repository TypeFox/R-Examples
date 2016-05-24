# For compiling the Fortran file:
# R CMD SHLIB vred.f


createRtopObject = function(observations, predictionLocations,
   formulaString, params=list(), ainfo, areas, overlapObs, overlapPredObs, 
   ...) {
  dots = list(...)
  if (inherits(observations,"rtop")) {
# Updating object with parameters
    object = observations
    object$params = getRtopParams(params, 
          observations = object$observations, formulaString = 
            ifelse(missing(formulaString), object$formulaString, formulaString), ...)
    return(object)
  }
  object = list()
  if (missing(observations) && !missing(ainfo) && !missing(areas)) {
    observations = SpatialPointsDataFrame(areas[ainfo$obsinc,],data=ainfo[ainfo$obsinc,])
    predictionLocations = SpatialPointsDataFrame(areas[!ainfo$obsinc,],data=ainfo[!ainfo$obsinc,])
  } else if ((!missing(observations) || !missing(predictionLocations)) && (!missing(ainfo) || !missing(areas))) {
    stop("ainfo and areas cannot be given as arguments together with observations and/or predictionLocations, please check documentation")
  } else {
    if (!inherits(observations,"SpatialPolygonsDataFrame") && "obs" %in% names(dots)) 
      if (!inherits(observations, "STS")) observations = SpatialPolygonsDataFrame(observations,data = dots$obs)
    if (!missing(predictionLocations) && !inherits(predictionLocations,"SpatialPolygonsDataFrame") && "pred" %in% names(dots)) 
      if (!inherits(predictionLocations, "STS"))predictionLocations = SpatialPolygonsDataFrame(predictionLocations,data = dots$pred)
  }
  if (missing(observations)) stop("Observations are missing")
#  if (missing(predictionLocations)) stop("predictionLocations are missing")
  if (!"area" %in% names(observations) && inherits(observations,"SpatialPolygons")) {
     observations$area = sapply(slot(observations, "polygons"), function(i) slot(i, "area"))
  } else if (inherits(observations, "STS") && !"area" %in% names(observations@sp)) {
    observations@sp$area = sapply(slot(observations@sp, "polygons"), function(i) slot(i, "area"))
  }
  object$observations = observations
  if (!missing(predictionLocations)) {
    if (!"area" %in% names(predictionLocations) && 
          inherits(predictionLocations,"SpatialPolygonsDataFrame")) {
      predictionLocations$area = sapply(slot(predictionLocations, "polygons"), function(i) slot(i, "area"))
    } else if (!"area" %in% names(predictionLocations) && 
          inherits(predictionLocations,"SpatialPolygons")) {
      areas = sapply(slot(predictionLocations, "polygons"), function(i) slot(i, "area"))
      predictionLocations = SpatialPolygonsDataFrame(predictionLocations, 
          data = data.frame(area = areas), match.ID = TRUE)  
#    } else if (!"length" %in% names(observations) && inherits(predictionLocations,"SpatialLines")) {
#       predictionLocations$length = SpatialLinesLengths(predictionLocations)
    } else if (inherits(predictionLocations, "STS") && !"area" %in% names(predictionLocations@sp)) {
      predictionLocations@sp$area = sapply(slot(predictionLocations@sp, "polygons"), function(i) slot(i, "area"))      
    } 
    object$predictionLocations = predictionLocations
  }
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
  object$formulaString = formulaString
#  depVar = formulaString[[2]] else depVar = "obs"
  object$params = getRtopParams(newPar = params,formulaString = formulaString,observations = observations)  
  if (length(dots) >0) object = modifyList(object,dots)  
  if (object$params$nugget) {
    if (!missing(overlapObs) && !is.null(overlapObs)) {
      object$overlapObs = overlapObs
    } else object$overlapObs = findOverlap(observations, debug.level = object$params$debug.level)
    if (!missing(overlapPredObs) && !is.null(overlapPredObs)) {
      object$overlapPredObs = overlapPredObs
    } else if (!missing(predictionLocations)) 
        object$overlapPredObs = findOverlap(observations,predictionLocations, debug.level = object$params$debug.level)
  }
  class(object) = "rtop" 
  object
}


getRtopParams = function(params = list(), newPar = list(), observations, formulaString, ...){
  
  
  oldPar = params
  dots = list(...) 
  if (inherits(oldPar,"intamapParams") || inherits(newPar,"intamapParams")) intPar = TRUE else intPar = FALSE
  oClass = class(oldPar)
  nClass = class(newPar)
  if (inherits(oldPar,"rtopParams")) {
    params = oldPar
    oldPar = list()
  } else if (inherits(newPar,"rtopParams")) {
    params = newPar
    newPar = list()
  } else {
    params = getRtopDefaultParams(...)
  }
  if (length(grep("geoDist", names(oldPar))) > 0 |
      length(grep("geoDist", names(newPar))) > 0 |
      length(grep("geoDist", names(dots))) > 0) stop("geoDist is not used anymore, please use gDist")
  
  params = modifyList(params, oldPar)
  params = modifyList(params, newPar)
  gDist = ifelse("gDist" %in% names(dots), dots$gDist,
            ifelse("gDist" %in% names(newPar), newPar$gDist,
            ifelse("gDist" %in% names(oldPar), oldPar$gDist,FALSE)))
  if (gDist) {
    params$gDistEst = TRUE
    params$gDistPred = TRUE
  }

  if (!missing(observations) & !("parInit" %in% names(params))) {
    if (missing(formulaString)) {
      if ("obs" %in% names(observations)) { 
        formulaString = "obs ~ 1" 
      } else if ("value" %in% names(observations)) {
        formulaString = "value ~ 1" 
      } else if (length(names(observations@data)) == 1) {
        formulaString = paste(names(observations@data),"~ 1")      
      } else stop("getRtopParams: formulaString is missing and cannot be found from data")
      warning(paste("getRtopParams: formulaString missing, using",formulaString))      
    }
    params$parInit = findParInit(formulaString,observations,params$model)
  } else if (!("parInit" %in% names(params))) {
    params$parInit = findParInitDefault(params$model)
  }
  params = modifyList(params,dots)

  if (intPar) class(params) = c("rtopParams","intamapParams") else class(params) = "rtopParams"
  params
}




getRtopDefaultParams = function(parInit,
   model="Ex1",
   nugget = FALSE,
   unc = TRUE,
   rresol = 100,  # Resolution real areas
   hresol = 5,    # Resolution in x-direction rectangles
#   logtrans = FALSE, # Logtransform data
   cloud = FALSE,  # work with cloud variogram
#   cutoff,        # cutoff distance in variogram - better to set in ... in call to function,
   amul = 2,		 # amul - defines the number of areal bins within one order of magnitude
   dmul = 3,		 # dmul - defines the number of distance bins within one order of magnitude
   fit.method = 9, # ils - Defines the type of Least Square method for fitting of variogram
                 #       1 - least squares difference  - err = yobs-ymod
                 #       2 - Weighted least squares difference according to Cressie (1985) - err2=n(yobs/ymod-1)^2
                 #       6 - No weights
                 #       7 - gstat fitting (Nj/hj^2)
                 #       8 - opposite of weighted least squares difference according to Cressie (1985) - err2=n*(ymod/yobs-1)^2
                 #       9 - Neutreal WLS-method - err = min(err2,err3)
   gDistEst = FALSE, # use ghosh distance
   gDistPred = FALSE,
   maxdist = Inf,
   nmax = 10,
   hstype = "regular", # Sampling type for hypothetical areas
#   rstype = ifelse(!missing(observations) && inherits(observations,"SpatialLines"),"regular","rtop"), 
                  # Sampling type for real areas
   rstype = "rtop",
   nclus = 1,
   cnAreas = 100,
   clusType = NULL,
   wlim = 1.5,
   wlimMethod = "all",
   cv = FALSE,
   debug.level = 1,
   observations,
   formulaString
   ){



#if (!missing(observations) & missing(cutoff)) {
#  x = coordinates(observations)[, 1]
#  y = coordinates(observations)[, 2]
#  cutoff = (0.35 * sqrt((max(x) - min(x))^2 + (max(y) - min(y))^2)/100)
#}
list(model = model, nugget = nugget, unc = unc, 
     rresol = rresol, hresol=hresol, rstype = rstype, hstype = hstype, 
#     logtrans = logtrans, 
     cloud = cloud, 
#     cutoff = cutoff, 
     amul = amul, dmul = dmul,
     fit.method = fit.method, gDistEst = gDistEst, gDistPred = gDistPred, 
     maxdist = maxdist, nmax = nmax, nclus = nclus, cnAreas = cnAreas, clusType = clusType,
     wlim = wlim, wlimMethod = wlimMethod, cv = cv,
     debug.level = debug.level)
}

###########################################
findParInitDefault = function(model) {

#  parameters are: sill, range, nugget, fractal, weibull par
  parInit = data.frame(parl = c(1e-06,   1e-02, 1.0e-01, 1e-5, 1e-01),
                       paru = c(5.0e+02, 1.0e7, 1.0e+07,  1.5,   1.7))
  parInit$par0 = 10**(0.5*(log10(parInit$paru)+log10(parInit$parl)))

  if (model %in% c("Exp","Sph", "Gau")) {
    parInit = parInit[1:3,]
  } else if (model == "Sp1"){
    parInit = parInit[1:4,]
  } else if (model == "Ex1"){
    parInit = parInit
  } else if (model == "Fra") {
    parInit[2,] = c(1e-6,2,0.01)
    parInit = parInit[1:3,]
  } else {
    stop(paste("model",model,"not implemented"))
  }
  parInit
}

#########################################
findParInit = function(formulaString,observations,model) {
  if (inherits(observations, "STS")) {
    ntime = dim(observations)[2]
    observations = observations[sample(1:ntime, 20),]
    vario = rtopVariogram(observations, formulaString = formulaString)
    aObs = sapply(slot(observations@sp, "polygons"), function(i) slot(i, "area"))
  } else {
    vario = variogram (formulaString, observations)
    aObs = sapply(slot(observations, "polygons"), function(i) slot(i, "area"))
  }
  parInit = data.frame(parl=c(1:5),paru=1,par0 = 1)
  parInit[1,1] = min(vario$gamma)/10
  parInit[1,2] = max(vario$gamma)*500 
  parInit[2,1] = sqrt(min(aObs))/4
  parInit[2,2] = max(vario$dist)*10
  minla = min(aObs)
  maxla = (max(aObs)^1.5)*max(vario$gamma)
  parInit[3,1] = min(vario$gamma)*minla/100
  parInit[3,2] = max(vario$gamma)*maxla
  parInit[4,1] = 1e-5
  parInit[4,2] = 1.5
  parInit[5,1] = 0.1
  parInit[5,2] = 1.7
  if (model == "Ex1") {
    parInit[4,2] = 1
    parInit[5,2] = 1  
  }
  

  parInit[,3] = sqrt(parInit[,1]*parInit[,2])
  if (model %in% c("Exp","Sph", "Gau")) {
    parInit = parInit[1:3,]
  } else if (model == "Sp1"){
    parInit = parInit[1:4,]
  } else if (model == "Ex1"){
    parInit = parInit
  } else if (model == "Fra") {
    parInit[2,] = c(1e-6,2,0.01)
    parInit = parInit[1:3,]
  } else {
    stop(paste("model",model,"not implemented"))
  }
  parInit
}

