rtopSim.rtop = function(object, varMatUpdate = FALSE, beta = NA, largeFirst = TRUE, params = list(), ...) {
  params = getRtopParams(object$params, newPar = params, ...)
  nmax = params$nmax
  cv = params$cv
  maxdist = params$maxdist
  wlim = params$wlim
  wlimMethod = params$wlimMethod
  debug.level = params$debug.level
  variogramModel = object$variogramModel
  if (is.null(variogramModel)) stop("Cannot do simulations without a variogram model")
  if (length(object$observations) > 0 &(is.null(object$varMatObs) | varMatUpdate))
    object = varMat(object, varMatUpdate, ...)
  if (is.null(object$varMatPred) || diff(dim(object$varMatPred)) != 0 | varMatUpdate) {
    if (is.null(object$dPred) & !(params$gDistPred & !is.null(object$gDistPred))) {
      object$dPred = rtopDisc(object$predictionLocations, params = params)
    }
    if (params$gDistPred & is.null(object$gDistPred)) {
      object$gDistPred = gDist(object$dPred, params = params)
      varMatPred = object$varMatPred = varMat(object$gDistPred, params = params, variogramModel = variogramModel)
    } else {
      varMatPred = object$varMatPred = varMat(object$dPred, params = params, variogramModel = variogramModel)
    }
  }
  varMatPredObs = object$varMatPredObs
  varMatObs = object$varMatObs
  varMatPred = object$varMatPred
  predictions = object$predictionLocations
  observations = object$observations
  nobs = length(observations)
  predictionLocations = object$predictionLocations
  predictions = predictionLocations
  if (!is(predictions, "SpatialPolygonsDataFrame")) {
    aPred = sapply(slot(predictions, "polygons"), function(i) slot(i, "area"))
    predictions = SpatialPolygonsDataFrame(predictions, data = data.frame(area = aPred ))
  } else if (!"area" %in% names(predictions)) {
    predictions$area = sapply(slot(predictions, "polygons"), function(i) slot(i, "area"))
  }
  singMat = FALSE
  varInv = NULL
  for (isim in 1:params$nsim) {
    predictions$sim = NA
    if (length(dim(observations)) > 0 && dim(observations)[1] > 0) {
      obsall = observations@data
      obs = obsall[,as.character(object$formulaString[[2]])]
      obscors = coordinates(observations)
      #      if (params$unc && "unc" %in% names(observations)) {
      #        unc0 = observations$unc
      #      } else unc0 = array(0,nobs)
      nobs0 = dim(observations)[1]
    } else {
      obsall = NULL
      obs = NULL
      obscors = NULL
      nobs0 = 0
    }
    vPred = varMatPred
    vObs = varMatObs
    vPredObs = varMatPredObs
    ips = sample(length(predictions), length(predictions))
    if (largeFirst) {
      il = which(ips == order(predictions$area, decreasing = TRUE)[1])
      tmp = ips[il] 
      ips[il] = ips[1]
      ips[1] = tmp
    }
    vpo = 1:length(predictions)
    if (interactive() & debug.level) {
      pb = txtProgressBar(1, length(ips), style = 3)
    }
    print(paste0(isim, ". simulation of ", length(ips), " areas"))
    for (ip in 1:length(ips)) {
      inew = ips[ip]
      in2 = which(vpo == inew)
      nobs = length(obs)
      if (interactive() & debug.level) setTxtProgressBar(pb, ip)
      newcor = coordinates(predictionLocations[inew,])
      unc0 = array(0,nobs)
      if (nobs == 0) {
        if (is.na(beta)) stop("No observations found, beta (expected mean) has to be given")
        c0 = varioEx(sqrt(bbArea(bbox(predictionLocations[in2,]))), variogramModel)
        inewvar = varMatPred[inew,inew]
        obs = rnorm(1, beta, c0-inewvar)
        vObs = matrix(inewvar, nrow = 1, ncol = 1)
        vPredObs = vPred[inew, -inew, drop = FALSE]
      } else {
        
        mdist = sqrt(diff(range(obscors[,1]))^2 + diff(range(obscors[,2]))^2)
        wlim0 = wlim
        while (TRUE) {
          wlim0 = wlim0/1.05
          ret <- try(rkrige(obsall, obs, obscors, newcor, vObs, vPredObs[,in2, drop = FALSE], nmax, inew, cv, 
                     unc0, mdist, maxdist, singMat, varInv, wlim0, debug.level, 
                     wlimMethod, simul = TRUE), silent = TRUE)
          if (is(ret, "try-error")) print(ip)
          if (wlim0 < 1.05 || (!is(ret, "try-error") && ret$pred[2] <= 0)) break
        }
        if (!is(ret, "try-error")) {
          nneigh = ret$nneigh
          lambda = ret$lambda
          neigh = ret$neigh
        
          pred = ret$pred
          newval = rnorm(1, pred[1], sqrt(pred[2]))
          obs = c(obs, newval)
        } else {
          obs = c(obs, NA)
        }
          vObs = rbind(vObs, vPredObs[,in2])
          vObs = cbind(vObs, c(vPredObs[,in2], 0))
          vPredObs = vPredObs[,-in2, drop = FALSE]
          vPredObs = rbind(vPredObs, vPred[in2, -in2])
      }
      obscors = rbind(obscors, newcor)        
      vPred = vPred[-in2, -in2, drop = FALSE]
      vpo = vpo[-in2]      
    }
    if (interactive() & debug.level) close(pb)
    predictions$sim = obs[(nobs0+1) : length(obs)]
    names(predictions)[dim(predictions)[2]] = paste0("sim", isim)
  }
  object$simulations = predictions
  object
}
 

rtopSim.default = function(object = NULL, predictionLocations, 
                           varMatObs, varMatPredObs, varMatPred, 
                           variogramModel, ...) {
object = createRtopObject(object, predictionLocations, ...)
if (!missing(varMatObs)) object$varMatObs = varMatObs
if (!missing(varMatPredObs)) object$varMatPredObs = varMatPredObs
if (!missing(varMatPred)) object$varMatPred = varMatPred
object$variogramModel = variogramModel
rtopSim(object, ...)$simulations
}
