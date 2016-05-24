spatialPredict.psgp = function(object,...) {

    dots = list(...)
    iparams = getIntamapParams(object$params, ...)
    # PSGP parameters are stored in the variogram model data frame
    # This is a hack to remain compatible with intamap. The first column
    # is discarded as it is for 
    # text and used to flag the values as incorrect (to ensure
    # this is not used as a variogram model). The other 9
    # columns contain the log parameters, which we put back into a vector
    params = object$variogramModel
    psgpLogParams = c(as.numeric(params[1,2:9]), as.numeric(params[2,2:9]))
    
    # vario=array()
    
    # variogram type
    # Gau - 1
    # Exp - 2
    #
    
    # vario[1] = as.integer(object$variogramModel$model[2])
    #if(object$variogramModel$model[2] == "Gau") vario[1]=1
    #if(object$variogramModel$model[2] == "Exp") vario[1]=2
    #vario[2]=object$variogramModel$range[2]
    #vario[3]=object$variogramModel$psill[2]
    #vario[4]=object$variogramModel$psill[1]
    #vario[5]=object$variogramModel$beta[1]
    
    
    rotated = FALSE
    if (object$params$doAnisotropy && object$anisPar$doRotation && all(as.character(object$formulaString[[3]])=="1")){
      objTemp = object
      object$observations = rotateAnisotropicData(object$observations, object$anisPar)
      object$predictionLocations = rotateAnisotropicData(object$predictionLocations, object$anisPar)
      rotated = TRUE
    }
    
    #if (require(astonGeostats)) {
    nPred = nrow(coordinates(object$predictionLocations))
    nsim = ifelse("nsim" %in% names(dots),dots$nsim,0) 
    if ("nclus" %in% names(object$params) && nsim == 0 && nPred >= 5000 ) 
      nclus = iparams$nclus else nclus = 1
    if (nclus > 1) {
      if (!suppressMessages(suppressWarnings(require(doParallel))))
  	    stop("nclus is > 1, but package doParallel is not available")    

      clus <- c(rep("localhost", nclus))
      cl <- makeCluster(clus, type = "SOCK")
      registerDoParallel(cl, nclus)
      clusterEvalQ(cl, library(psgp))
      formulaString = object$formulaString
      observations = object$observations
      predictionLocations = object$predictionLocations
      variogramModel = object$variogramModel
#      clusterExport(cl, list("formulaString", "observations", "predictionLocations",
#           "variogramModel", "nmax", "nsim", "debug.level"))
     # split prediction locations:
      splt = sample(1:nclus, nPred, replace = TRUE)
      splt = rep(1:nclus, each = ceiling(nPred/nclus), length.out = nPred)
      newdlst = lapply(as.list(1:nclus), function(w) predictionLocations[splt == w,])
      nobject = object
      i <- 1
      pred <- foreach(i = 1:nclus) %dopar% {
        nobject$predictionLocations = newdlst[[i]]
        makePrediction(nobject, psgpLogParams)
      }
      stopCluster(cl)
      var1.pred = unlist(lapply(pred,FUN = function(pp) pp[[1]]))
      var1.var = unlist(lapply(pred,FUN = function(pp) pp[[2]]))
      object$predictions = SpatialPointsDataFrame(object$predictionLocations,
        data = data.frame(var1.pred = var1.pred,var1.var=var1.var))
    } else {  
      pred = makePrediction(object, psgpLogParams)
      object$predictions = SpatialPointsDataFrame(object$predictionLocations,
        data = data.frame(var1.pred = unlist(pred[1]),var1.var=unlist(pred[2])))
    }
    if (nsim > 0) {
      nmax = object$params$nmax
      object$predictions = cbind(object$predictions,krige(object$formulaString,object$observations, 
             object$predictionLocations,object$variogramModel,nsim=nsim,nmax = nmax,debug.level = object$params$debug.level))
    }
    if (rotated) {
      object$observations = objTemp$observations
      object$predictionLocations = objTemp$predictionLocations
      object$predictions@coords = coordinates(object$predictionLocations)
      object$predictions@bbox = bbox(object$predictionLocations)
      proj4string(object$predictions) = proj4string(object$predictionLocations)
    }
    names(object$predictions) = c("var1.pred","var1.var")
    object
}

