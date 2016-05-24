
# varMat.rtop hardly uses the functions further below, should be shortened

varMat.rtop = function(object, varMatUpdate = FALSE, params = list(), ...) {
  params = getRtopParams(object$params,  newPar = params, ...)
  observations = object$observations
  if (is(observations, "STSDF")) observations = observations@sp
  nObs = dim(observations@data)[1]
  predictionLocations = object$predictionLocations
  if (is(predictionLocations, "STSDF")) predictionLocations = predictionLocations@sp
  variogramModel = object$variogramModel
  lgDistPred = params$gDistPred
  maxDist = params$maxDist
  debug.level = params$debug.level
  aObs = sapply(slot(observations, "polygons"), function(i) slot(i, "area"))

  obsComp = FALSE
  predComp = FALSE
  if ("varMatObs" %in% names(object) && !varMatUpdate) 
    if (!identical(attr(object$varMatObs, "variogramModel"),variogramModel) | 
        !nObs == dim(object$varMatObs)[1] |  
         ("varMatPredObs" %in% names(object) && (
            dim(object$varMatPredObs)[2] != dim(predictionLocations)[1] |
            !identical(attr(object$varMatObs, "variogramModel"),variogramModel))))
               varMatUpdate = TRUE
  if (params$cv && "varMatObs" %in% names(object) && !varMatUpdate) return(object)      
  if (!"varMatObs" %in% names(object) | varMatUpdate) {
    if (!"dObs" %in% names(object) && !(lgDistPred && "gDistObs" %in% names(object))) 
        object$dObs = rtopDisc(observations, params = params) 
    if (lgDistPred) {
      if ("gDistObs" %in% names(object)) {
        gDistObs = object$gDistObs
      } else {
        dObs = object$dObs
        object$gDistObs = gDistObs = gDist(dObs, params = params)
      }
      if (!is.null(params$nclus) && params$nclus > 1 && nObs > params$cnAreas && requireNamespace("parallel")) {
        cl = rtopCluster(params$nclus, type = params$clusType)
        varMatObs = matrix(unlist(parallel::parLapply(cl, 1:nObs, fun = function(x, gDistObs, variogramModel) 
          mapply(gDistObs[, x],
                 FUN = function(y) varioEx(y, variogramModel)), 
                           gDistObs = gDistObs, variogramModel = variogramModel)), 
                           nrow = nObs, ncol = nObs)
      } else {
        varMatObs = matrix(mapply(gDistObs, FUN = varioEx, 
                          MoreArgs=list(variogramModel)),
                          nrow = nObs,ncol = nObs)
      }
      vDiagObs = diag(varMatObs)
      for (ia in 1:(nObs-1)) {
        for (ib in (ia+1):nObs) {
          varMatObs[ia,ib] = varMatObs[ia,ib] - 0.5*(vDiagObs[ia] + vDiagObs[ib])
          varMatObs[ib,ia] = varMatObs[ia,ib]
        }
      }
      object$varMatObs = varMatObs
    } else {
      if ("dObs" %in% names(object)) dObs = object$dObs
      object$varMatObs = varMat(dObs,coor1 = coordinates(observations), 
          variogramModel = variogramModel, debug.level = debug.level, newPar = params)
    }
    attr(object$varMatObs, "variogramModel") = variogramModel
    obsComp = TRUE
  }
  if (!params$cv && "predictionLocations" %in% names(object) && 
        (!"varMatPredObs" %in% names(object) | varMatUpdate)) {
#    ftype = ifelse(inherits(predictionLocations,"SpatialPolygons"),"polygons","lines")
    varMatObs = object$varMatObs
    vDiagObs = diag(varMatObs)
    ftype = "polygons"
    nPred = length(sapply(slot(predictionLocations, ftype), function(i) slot(i, "ID")))
    if (!"dPred" %in% names(object) && !(lgDistPred && "gDistPred" %in% names(object))) 
      object$dPred = rtopDisc(predictionLocations, params = params) 
    if (!lgDistPred || !all("gDistPredObs" %in% names(object) &&  "gDistPred" %in% names(object))) {
      dObs = object$dObs
      dPred = object$dPred
    }
    if (lgDistPred) {
      if ("gDistPred" %in% names(object)) {
        gDistPred = object$gDistPred
      } else object$gDistPred = gDistPred = gDist(dPred, diag=TRUE, params = params)
      if ("gDistPredObs" %in% names(object)) {
        gDistPredObs = object$gDistPredObs
      } else object$gDistPredObs = gDistPredObs = gDist(dObs,dPred, params = params)

      print("Creating prediction semivariance matrix. This can take some time.")
      object$varMatPred = varMatPred = 
                matrix(mapply(FUN = varioEx,gDistPred,MoreArgs=list(variogramModel)),
                          nrow = nPred,ncol = 1)
      
      if (!is.null(params$nclus) && params$nclus > 1 && nObs > params$cnAreas && requireNamespace("parallel")) {
        cl = rtopCluster(nclus = params$nclus, type = params$clusType)

        varMatPredObs = matrix(unlist(parallel::parLapply(cl, 1:nPred, fun = function(x, gDistPredObs, variogramModel) 
          mapply(gDistPredObs[, x],
                 FUN = function(y) varioEx(y, variogramModel)), 
          gDistPredObs = gDistPredObs, variogramModel = variogramModel)), 
          nrow = nObs, ncol = nPred)      
      } else {
        varMatPredObs = matrix(mapply(FUN = varioEx,gDistPredObs,MoreArgs=list(variogramModel)),
                          nrow = nObs,ncol = nPred)
      }
      if (is.null(dim(varMatPred)) || dim(varMatPred)[1] != dim(varMatPred)[2]) 
          vDiagPred = varMatPred else vDiagPred = diag(varMatPred)
      for (ia in 1:nObs) {
        for (ib in 1:nPred) 
           varMatPredObs[ia,ib] = varMatPredObs[ia,ib] - 0.5*(vDiagObs[ia] + vDiagPred[ib])
      }
      object$varMatPredObs = varMatPredObs    

    } else {
   # Do full integration over variograms
      object$varMatPred = varMat(dPred,coor1 = coordinates(predictionLocations), 
                    diag=TRUE, variogramModel = variogramModel, debug.level = debug.level, newPar = params)
      object$varMatPredObs = varMat(dObs,dPred,
          coor1 = coordinates(observations),coor2 = coordinates(predictionLocations),
          variogramModel = variogramModel, sub1 = diag(object$varMatObs), sub2 = object$varMatPred,
          debug.level = debug.level, newPar = params)
    }
    predComp = TRUE
  }
  if (params$nugget) {
    if (obsComp) {
      if ("overlapObs" %in% names(object)) {
        overlapObs = object$overlapObs
      } else object$overlapObs = overlapObs = findOverlap(observations,observations, params = params)
      fObs = matrix(rep(aObs,nObs),ncol=nObs)
      sObs = t(fObs)
      nuggObs = matrix(mapply(FUN = nuggEx,
            (1/fObs + 1/sObs -2*overlapObs/(fObs*sObs))/2,
             MoreArgs = list(variogramModel = variogramModel)),ncol = nObs)
      diag(nuggObs) = 0
      object$varMatObs = object$varMatObs + nuggObs
    }
    if (predComp) {
      if ("overlapPredObs" %in% names(object)) {
        overlapPredObs = object$overlapPredObs
      } else object$overlapPredObs = overlapPredObs = findOverlap(observations,predictionLocations, params = params)
      aPred = sapply(slot(predictionLocations, "polygons"), function(i) slot(i, "area"))
      fPredObs = matrix(rep(aObs,nPred),ncol=nPred)
      sPredObs = t(matrix(rep(aPred,nObs),ncol = nObs))
      nuggPredObs = matrix(mapply(FUN = nuggEx,
            (1/fPredObs + 1/sPredObs -2*overlapPredObs/(fPredObs*sPredObs))/2,
            MoreArgs = list(variogramModel = variogramModel)),ncol= nPred)
      object$varMatPredObs = object$varMatPredObs + nuggPredObs
    }
  }
  if ("varMatPredObs" %in% names(object)) attr(object$varMatPredObs, "variogramModel") = variogramModel
  if ("varMatPred" %in% names(object)) attr(object$varMatPred, "variogramModel") = variogramModel
  object
}
    

varMat.matrix = function(object, variogramModel, diag = FALSE, sub1, sub2, ...) {
  ndim = dim(object)[1] 
  mdim = dim(object)[2]
  varMatrix = matrix(mapply(FUN = varioEx,object,MoreArgs=list(variogramModel)),
                          nrow = ndim,ncol = mdim)
  if (diag) sub1 = sub2 = diag(varMatrix)
  if (!missing(sub1) & !missing(sub2)) {  
    for (ia in 1:ndim) {
      for (ib in 1:mdim) if (!(diag & ia == ib)) varMatrix[ia,ib] = varMatrix[ia,ib] - 0.5*(sub1[ia] + sub2[ib])
    }
  }     
  varMatrix
  
}


varMat.SpatialPolygonsDataFrame = function(object,object2 = NULL,...) {
  if (is(object2,"SpatialPolygonsDataFrame")) object2 = as(object2,"SpatialPolygons")
  varMat(as(object,"SpatialPolygons"), object2, ...)

}




varMat.SpatialPolygons = function(object,object2 = NULL,variogramModel,
       overlapObs, overlapPredObs, ...) {
  varMatDefault(object,object2,variogramModel,
     overlapObs, overlapPredObs, ...) 
}


varMatDefault = function(object1,object2 = NULL,variogramModel,
     overlapObs, overlapPredObs,  ...) {
  params = getRtopParams(...)
  d1 = rtopDisc(object1, params) 
  if (!is.null(object2)) d2 = rtopDisc(object2, params) 
  if (params$gDistPred) {
    gDist1 = gDist(d1, params = params)
# calling varMat.matrix
    varMatObs = varMat(gDist1, variogramModel = variogramModel, params = params, diag = TRUE, ...)
  } else {
    varMatObs = varMat(d1, variogramModel = variogramModel, params = params, ...)
  }  

  if (is.null(object2)) return(varMatObs)
  
  if (params$gDistPred & !is.null(object2)) {
    gDistPred = gDist(d2, diag = TRUE, params = params)
# Calling varMat.matrix
    varMatPred = varMat(gDistPred, params = params, variogramModel = variogramModel, ...)
    gDistPredObs = gDist(d1, d2, params = params)
    varMatPredObs = varMat(gDistPredObs,sub1 = diag(varMatObs),sub2 = varMatPred, params = params, 
                           variogramModel = variogramModel, ...)
  } else {
# Calling varMat.list
    varMatPred = varMat(d2,diag=TRUE, params = params, variogramModel = variogramModel, )
    varMatPredObs = varMat(d1,d2,sub1 = diag(varMatObs),sub2 = varMatPred, params = params, 
                           variogramModel = variogramModel, ...)
  }
  if (params$nugget) {
    if (missing(overlapObs)) 
      overlapObs = findOverlap(object1,object1, params = params)
    if (missing(overlapPredObs))
      overlapPredObs = findOverlap(object1,object2, params = params)
    
    aObs = sapply(slot(object1, "polygons"), function(i) slot(i, "area"))
    aPred = sapply(slot(object2, "polygons"), function(i) slot(i, "area"))
    
    nObs = length(aObs)
    nPred = length(aPred)
    
    fObs = matrix(rep(aObs,nObs),ncol=nObs)
    sObs = t(fObs)
    fPredObs = matrix(rep(aPred,nPred),ncol=nPred)
    sPredObs = t(matrix(rep(aPred,nObs),ncol = nObs))
    nuggObs = matrix(mapply(FUN = nuggEx,
            (1/fObs + 1/sObs -2*overlapObs/(fObs*sObs))/2,
             MoreArgs = list(variogramModel = variogramModel)),ncol = nObs)
    nuggPredObs = matrix(mapply(FUN = nuggEx,
            (1/fPredObs + 1/sPredObs -2*overlapPredObs/(fPredObs*sPredObs))/2,
            MoreArgs = list(variogramModel = variogramModel)),ncol = nPred)
    object$varMatObs = object$varMatObs - nuggObs
    object$varMatPredObs = object$varMatPredObs - nuggPredObs
  }
  attr(varMatObs, "variogramModel") = variogramModel
  attr(varMatPredObs, "variogramModel") = variogramModel
  attr(varMatPred, "variogramModel") = variogramModel
  return(list(varMatObs = varMatObs,varMatPred = varMatPred,varMatPredObs = varMatPredObs))  
}


# object and object2 (as lists) are here discretized areas
# coor1 and coor2 are coordinates of the areas, used for maximum distance
varMat.list = function(object, object2=NULL, coor1, coor2, maxdist = Inf, 
              variogramModel, diag=FALSE, sub1, sub2, debug.level = 1, ...) {
  params = getRtopParams(...)
  d1 = object
  d2 = object2
  if (is.null(d2)) {
    equal = TRUE
    d2 = d1
    if (!missing(coor1)) coor2 = coor1
  } else equal = FALSE
  if (diag) {
    lmat = mapply(vred,a2 = d1,a1 = d1,MoreArgs = list(vredTyp="ind",variogramModel = variogramModel))
    return(lmat)
  }

  ndim = length(d1)
  mdim = length(d2)
  varMatrix = matrix(-999,nrow = ndim,ncol = mdim)

  if (!is.null(params$nclus) && params$nclus > 1 && length(d1) + length(d2) > params$cnAreas && requireNamespace("parallel")) {
#    cl = rtopCluster(params$nclus, {require(sp); vred = rtop:::vred}, type = params$clusType)
    cl = rtopCluster(params$nclus, type = params$clusType)
    if (missing(coor1)) coor1 = NULL
    if (missing(coor2)) coor2 = NULL

    
    fun = function(ia, d1, d2, coor1, coor2, equal, maxdist) {    
      ndim = length(d1)
      mdim = length(d2)
      a1 = coordinates(d1[[ia]])
      ip1 = dim(a1)[1]
      first = ifelse(equal,ia,1)
      lorder = c(first:mdim)
      if (!is.null(coor1) && ! is.null(coor2) && maxdist < Inf) 
          lorder = lorder[spDistsN1(coor2[first:mdim,],coor1[ia,]) < maxdist]
        if (length(lorder) > 0) {
          ld = d2[lorder]
          lmat = mapply(vred,a2 = ld,MoreArgs = list(vredTyp="ind",a1 = a1,variogramModel = variogramModel))
        } else lmat = -999  
      list(lmat, lorder)
    }    
    vmll = parallel::clusterApply(cl, 1:length(d1), d1 = d1, d2 = d2, coor1 = coor1, coor2 = coor2, 
                       equal = equal, maxdist = maxdist, fun = fun)

    for (ia in 1:length(vmll)) {
      lmat = vmll[[ia]][[1]]
      lorder = vmll[[ia]][[2]]
      if (!equal) varMatrix[ia,] = lmat else {
        varMatrix[ia,lorder] = lmat
        varMatrix[lorder,ia] = lmat
      }
    }
  } else {
    t0 = proc.time()[[3]]
    for (ia in 1:ndim) {
      t1 = proc.time()[[3]]
      a1 = coordinates(d1[[ia]])
      ip1 = dim(a1)[1]
      first = ifelse(equal,ia,1)
      lorder = c(first:mdim)
      if (!missing(coor1) && ! missing(coor2) && maxdist < Inf) 
           lorder = lorder[spDistsN1(coor2[first:mdim,],coor1[ia,]) < maxdist]
      if (length(lorder) > 0) {
        ld = d2[lorder]
        lmat = mapply(vred,a2 = ld,MoreArgs = list(vredTyp="ind",a1 = a1,variogramModel = variogramModel))
#        lmat = lvmat[[1]]
        if (!equal) varMatrix[ia,] = lmat else {
          varMatrix[ia,lorder] = lmat
          varMatrix[lorder,ia] = lmat
        }
      }
#        print(acdf)
      t2 = proc.time()[[3]]
      if (debug.level > 0) print(paste("varMat - Finished element ",ia," out of ", ndim," in ", round(t2-t1,3),
            "seconds - totally", round(t2-t0)," seconds"))
    }
  }
    if (equal & variogramModel$model != "Gho") {
      vDiag = diag(varMatrix)
      for (ia in 1:(ndim-1)) {
        for (ib in (ia+1):ndim) {
          varMatrix[ia,ib] = varMatrix[ia,ib] - 0.5*(vDiag[ia] + vDiag[ib])
          varMatrix[ib,ia] = varMatrix[ia,ib]
        }
    }
  } else if (!missing(sub1) & !missing(sub2)) {  
    for (ia in 1:ndim) {
      for (ib in 1:mdim) varMatrix[ia,ib] = varMatrix[ia,ib] - 0.5*(sub1[ia] + sub2[ib])
    }
  }     
  attr(varMatrix, "variogramModel") = variogramModel
  varMatrix
}
  


varMat.STS = function(object, object2 = NULL, variogramModel, overlapObs, overlapPredObs, ...) {
  if (is(object2,"STS")) object2 = object2@sp
  object = object@sp
  if (!is(object, "SpatialPolygons") || (!is.null(object2) && !inherits(object2, "SpatialPolygons")))
    stop(paste("Cannot create covariance matrix from objects of class", class(object), class(object2))) 
  varMat(object, object2, variogramModel, overlapObs, overlapPredObs, ...)
}