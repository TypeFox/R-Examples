  
blockPredict = function(object,...) {
  object = spatialPredict.block(object,...)
  return(object)
} 
 

spatialPredict.block = function(object,...) {
#spablock = function(object,...) {
# Define which methods can actually handle block kriging directly
  blockMethods = c("automap","idw")
# Extract some of the parameters and variables from the object
# - for easier access
  params = getIntamapParams(object$params, ...)
  blockWhat = object$blockWhat
  dots = list(...)
  block = params$block
  nmax = params$nmax
  observations = object$observations
  formul = object$formulaString
  predictionLocations = object$predictionLocations
  outputWhat = object$outputWhat
  if (!inherits(predictionLocations,"SpatialPolygons") &
      !inherits(predictionLocations,"SpatialLines") &
      !inherits(predictionLocations,"SpatialPixels") &
      !inherits(predictionLocations,"SpatialGrid") &
      !length(block) >0 ) {
      warning("prediction locations are not grid, blocks, polygons or lines")
      warning("calling point prediction")
      return(spatialPredict(object,...))
  }
  
  if (!("variogramModel" %in% names(object)) &
      !inherits(object,"idw") & !inherits(object,"copula")) object = estimateParameters(object, ...)
  if (params$processType %in% c("gaussian","logNormal") &
      sum(class(object) %in% blockMethods) >= 1) {
# Check if we might be able to use an analytical solution - 
# Possible if we want a prediction for Gaussian or logNormal distribution, only 1 cluster
# Only possible for certain methods - now hard coded into blockMethods
    if (inherits(object,"idw")) return(spatialPredict(object,...))
    form = object$formulastring
    depVar = as.character(form[2])
# We have to discuss whether it makes sense to have the test of logNormal here, 
# as it can easily create some inconsistencies.
# We have to be sure that the variogram is also found from the logarithmised data
    if (params$processType == "logNormal") observations[[depVar]] = log(observations[[depVar]])
      predictions = krige(object$formulaString,observations, predictionLocations,object$variogramModel,block = block,nmax = nmax)
    if (params$processType == "logNormal") {
      predictions$var1.pred = exp(predictions$var1.pred + predictions$var1.var/2) 
      warning("Note that the back-calculated lognormal predictor is only approximate, as the Lagrange parameter is not included")
    }
    print("performed ordinary block kriging")
    object$predictions = predictions
  }
  
# Point predictions/simulations are necessary
    if ("MOK" %in% names(object$outputWhat) | "IWQSEL" %in% names(object$outputWhat)
      | inherits(object,"yamamoto") | length(names(blockWhat)) > 0
      | !inherits(object,blockMethods)) {
# Create a new object for point simulation
    pointObject = object
#   Create a grid
    if (!inherits(predictionLocations, "SpatialPolygons")) {
      if (length(block) > 0) {
        lbox = vector("list", length(predictionLocations))      
        coords = coordinates(predictionLocations)
        for (ii in 1:length(predictionLocations)) {
          xy = coords[ii,]
          xmin = xy[1]-block[1]/2
          xmax = xy[1]+block[1]/2
          ymin = xy[2]-block[ifelse(length(block) == 2, 2, 1)]/2
          ymax = xy[2]+block[ifelse(length(block) == 2, 2, 1)]/2
          Sr1 = Polygon(cbind(c(xmin,xmax,xmax,xmin,xmin),c(ymin,ymin,ymax,ymax,ymin)))
          lbox[[ii]] = Polygons(list(Sr1), paste("i",ii,sep = ""))
        }
        SpP = SpatialPolygons(lbox, 1:length(predictionLocations), 
              proj4string = CRS(proj4string(predictionLocations)))      
      } else if (inherits(predictionLocations, "SpatialPixels") | 
                 inherits(predictionLocations, "SpatialPixels")) {
        SpP = as(predictionLocations, "SpatialPolygons")      
      } else stop("predictionLocations are not SpatialPolygons or gridded and block size not given, block prediction not possible ")
    } else SpP = predictionLocations

    pointObject$predictionLocations = findGrid(predictionLocations = SpP,
          params = params)

# Do we need point predictions or simulations?
# Prediction only if only need mean
    
    nsim = ifelse(all(names(outputWhat)=="mean") & blockWhat == "none",0, params$nsim) 
    vmod = object$variogramModel
    nmax = ifelse(params$nmax == Inf & dim(coordinates(pointObject$predictionLocations))[1] > 200,20,params$nmax)
    pointObject = spatialPredict(pointObject,nsim=params$nsim,nmax = nmax, ...)
    pointObject$predictions = pointObject$predictions[,grep("sim",names(pointObject$predictions))]
    object$pointPredictions = pointObject$predictions
    object$pointLocations = pointObject$predictionLocations
    object = spatialAggregate(object, SpP)
    if (object$params$debug.level < 2) object$pointPredictions = 
                  "pointPredictions deleted from object, debug.level < 2"
  }
  return(object)

}


  
spatialAggregate = function(object, SpP) {
  predictionLocations = object$predictionLocations
  pointLocations = object$pointLocations  
  if ("predictions" %in% names(object)) {
    predictions = object$predictions
  } else {
    predictions = predictionLocations
  }
    
  pointPredictions = object$pointPredictions
  coor = SpatialPoints(pointPredictions)
  sims = pointPredictions[,grep("sim",names(pointPredictions))>0]
  
  params = object$params
  thresh = params$thresh
  outputWhat = object$outputWhat
  blockWhat = object$blockWhat
  if ("mean" %in% names(outputWhat)) {
     if (blockWhat == "none") blockWhat = list(mean = TRUE) else blockWhat$mean = TRUE
  }
  predAggr = aggregate(sims,by=SpP,mean)
  if ("data" %in% names(getSlots(class(predictions)))) {
    predictions@data = data.frame(predictions@data,predAggr@data)
  } else predictions = SpatialDataFrame(predictions,predAggr@data)
  if (length(blockWhat) > 0 && blockWhat != "none") {
    for (ib in 1:length(blockWhat)) {
      what = blockWhat[ib]
      if (names(what) == "fat") {
        thresh = what[[1]]
        fatf = function(arr,thresh) sum(I(arr>thresh))/length(arr)  
        fatx = aggregate(sims,SpP,FUN = fatf,thresh=thresh)@data
        vmean = rowMeans(fatx[,-1])
        vvar = apply(fatx[,-1],MARGIN=1,FUN=function(arr) var(arr))
        vname = paste("fat",what[[1]],sep="")
        vnamevar = paste("fatVar",what[[1]],sep="")
      } else if (names(what) == "blockMax" && what[[1]]) {
        bmax = aggregate(sims,by=SpP,FUN = max)@data
        vmean = rowMeans(bmax[,-1])
        vvar = apply(bmax[,-1],MARGIN=1,FUN=function(arr) var(arr))
        vname = "blockMax"
        vnamevar = "blockMaxVar"
      } else if (names(what) == "blockMin" && what[[1]]) {
        bmin = aggregate(sims,by=SpP,FUN = min)@data
        vmean = rowMeans(bmin[,-1])
        vvar = apply(bmin[,-1],MARGIN=1,FUN=function(arr) var(arr))
        vname = "blockMin"
        vnamevar = "blockMinVar"
      } else if (names(what) == "mean" && what[[1]]) {
        bmean = aggregate(sims, by = SpP, FUN = mean)@data
        vmean = rowMeans(bmean[,-1])
        vvar = apply(bmean[,-1],MARGIN=1,FUN=function(arr) var(arr))
        vname = "blockMean"
        vnamevar = "blockMeanVar"
      }
      predictions@data[vname] = vmean        
      predictions@data[vnamevar] = vvar             
    }
  }
  object$predictions = predictions
  object
}  

  
                  
findGrid = function(predictionLocations, params) {
#  newdata needs to be SpatialPointsDataFrame.
#  3 options:
#  1) If sampleSubregions is false then sample the whole region of newdata
#  2) If sampleSubregions is false subregions denominated the string "regCode" 
#     are sampled
#  3) If block is given, 
#  If both n and cellsize is given, cellsize is used, unless the sampling
#  will be done in subregions, where cellsize does not make sense
#  If cellmin gives empty polygons, sMin gives a minimum number of samples from
#  each polygon   
sMin = params$sMin
ncell = params$ngrid
cellsize = params$cellsize
sampleSubregions = params$subSamp

  if (!is.null(sampleSubregions) && sampleSubregions) {
    if (!missing(cellsize)) ncell = bbArea(bbox(predictionLocations))/(cellsize*cellsize)
    if (ncell < sMin) ncell = sMin
    ids = sapply(slot(predictionLocations, "polygons"), function(i) slot(i, "ID"))
    for (i in 1:length(predictionLocations@polygons)) {
      ldata = predictionLocations@polygons[i][[1]]
      nl = length(ldata@Polygons)
      cArea = 0
    # Summing up areas of all subpolygons in a polygon
      for (ip in 1:nl) cArea = cArea + ldata@Polygons[[ip]]@area
      aR = spsample(ldata,type="regular", n=ncell)
      naR = length(coordinates(aR)[,1])      
      if (is(predictionLocations,"SpatialPolygonsDataFrame")) {
        idl = ids[i]
        idNew = as.character(rep(idl,naR))
      } else idNew = as.character(rep(i,naR))
      predGridNew = SpatialPointsDataFrame(aR,data=data.frame(id = idNew))
# 
      if (i == 1) {
        id = idNew
        predGrid = predGridNew
      } else {
        id = c(id,idNew)
        predGridNew = SpatialPointsDataFrame(aR,data=data.frame(id = idNew))                     
        predGrid = rbind(predGrid,predGridNew)
      }
    }
  } else {  
    if (!is.null(cellsize)) {
      predLoc = spsample(predictionLocations,type="regular",cellsize=cellsize) 
    } else {
      predGrid = spsample(predictionLocations,ncell*length(predictionLocations),type="regular")
    }
  }

if (!is.na(proj4string(predictionLocations))) proj4string(predGrid) = proj4string(predictionLocations) 
predGrid
}

