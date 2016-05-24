
rtopFitVariogram.rtop = function(object, params = list(), ...) {

  params = getRtopParams(object$params, newPar = params, ...)
  if (params$cloud && !"variogramCloud" %in% names(object) ||
     !params$cloud && !"variogram" %in% names(object))
    object = rtopVariogram(object,...)

  if (params$nugget) {
    if (!"overlapObs" %in% names(object)) 
      object$overlapObs = overlapObs = findOverlap(observations,observations)
  }  
  if (params$cloud) {
    vario = object$variogramCloud
    observations = object$observations
    if (params$gDistEst && "gDistObs" %in% names(object)) {
      dists = object$gDistObs
    } else if ("dObs" %in% names(object)) {
      dists = object$dObs
    } else dists = NULL
    if (params$nugget) aOver = object$overlapObs else aOver = NULL
  } else {
    vario = object$variogram
    observations = NULL
    if (params$gDistEst && "gDistBin" %in% names(object)) {
      dists = object$gDistBin
    } else if ("dBin" %in% names(object)){
      dists = object$dBin
    } else dists = NULL
    if (params$nugget) aOver = findVarioOverlap(vario) else aOver = NULL
  }
  varioFit = rtopFitVariogram(object = vario, observations = observations,
             dists = dists, aOver = aOver, params = params, mr = TRUE, ...)
  object$variogramModel = varioFit$variogramModel
  object$varFit = varioFit$varFit
  varioFit = varioFit[-which(names(varioFit) == "variogramModel")]
  varioFit = varioFit[-which(names(varioFit) == "varFit")]
# Adding gdist or discretized areas if created
  object = modifyList(object, varioFit)
  object
}

# Dists can be
# dObs
# dsBin

rtopFitVariogram.rtopVariogram = function(object, observations, dists = NULL, params=list(), mr = FALSE, aOver = NULL, ...){
  vario = object
  obj = list()
  if (!inherits(params,"rtopParams")) params = getRtopParams(params, ...)
  obj$params = params
  if (is.null(dists)) {
    obj$dBin = dists = rtopDisc(object, params = params)
  }
  if (params$gDistEst && !is.matrix(dists)) obj$gDistBin = dists = gDist(dists, params = params)
  if (params$nugget & is.null(aOver)) aOver = findVarioOverlap(vario) 

  if (params$model == "Ex1") {
    implicit = function(pars) (2*pars[4] + pars[5]) > 1
  } else implicit = NULL
  scres = sceua(objfunc,params$parInit[,3],params$parInit[,1],params$parInit[,2],varioIn = object,
       dists = dists, aOver = aOver, gDistEst = params$gDistEst,model = params$model,resol = params$hresol,
       fit.method = params$fit.method, implicit = implicit, ...)
  bestPar = scres$par
  fit = scres$value
  vf = objfunc(bestPar,varioIn = vario, dists = dists, aOver = aOver,
                   gDistEst = params$gDistEst,last = TRUE,model = params$model,resol = params$hresol,...)
  varFit = vf$varFit
  errSum = vf$errSum

  variogramModel = list(model = params$model,params = bestPar)
  class(variogramModel) = "rtopVariogramModel"
  attr(variogramModel,"SSErr") = errSum
  attr(variogramModel,"criterion") = fit
  if (!params$nugget) variogramModel$params[3] = 0
  if (mr) {
# Here adding to the object created further up
    obj$variogramModel = variogramModel
    obj$varFit = varFit
    obj
  } else list(variogramModel = variogramModel,varFit = varFit)
}

# Dists can be
# gDists
# discAreas

rtopFitVariogram.rtopVariogramCloud = function(object, observations, dists = NULL, aOver = NULL, params=list(), mr = FALSE, ...) {
  vario = object
  obj = list()
  if (!inherits(params,"rtopParams")) params = getRtopParams(params, ...)
  obj$params = params
  if (is.null(dists))
    obj$dObs = dists = rtopDisc(observations,params = params,...)
  if (params$gDistEst && is.list(dists))
    obj$gDistObs = dists = gDist(dists, params = params, ...)

  if (params$model == "Ex1") {
    implicit = function(pars) (2*pars[4] + pars[5]) > 1
  } else implicit = NULL
  scres = sceua(objfunc, params$parInit[,3], params$parInit[,1], params$parInit[,2], varioIn = vario,
         dists = dists, aOver = aOver, gDist = params$gDistEst,model = params$model, 
         fit.method = params$fit.method, debug.level = params$debug.level, ...)
  bestPar = scres$par
  vf = objfunc(bestPar,varioIn = vario, dists = dists, aOver = aOver,
               gDistEst = params$gDistEst, last = TRUE, model = params$model,
               debug.level = params$debug.level, ...)
  varFit = vf$varFit
  errSum = vf$errSum
  variogramModel = list(model = params$model, params = bestPar)
  class(variogramModel) = "rtopVariogramModel"
  attr(variogramModel, "SSErr") = errSum
  if (mr) {
    obj$variogramModel = variogramModel
    obj$varFit = varFit
    obj
  } else list(variogramModel ,varFit = varFit)
}


rtopFitVariogram.SpatialPointsDataFrame = function(object, params=list(),...) {
  if (!inherits(params,"rtopParams")) params = getRtopParams(params, ...)
  vario = rtopVariogram(object,params,...)
  rtopFitVariogram(vario = vario, observations = object, params = params,...)
}

rtopFitVariogram.SpatialPolygonsDataFrame = function(object, params=list(),...) {
  if (!inherits(params,"rtopParams")) params = getRtopParams(params, ...)
  vario = rtopVariogram(object,params = params,...)
  rtopFitVariogram(vario = vario, observations = object, params = params,...)
}






objfunc = function(pars, varioIn, gDistEst=FALSE, dists, aOver = NULL, model="Ex1",
                  bu, bl, fit.method = 8, debug.level = 0, resol = 5, nd =100,
                  last = FALSE, ...) {
# Debug = 0 means no output
# Debug = 1 means output for every nd iteration
# Debug = 2 means output for every iteration
# Debug = 3 means output for every element in variogram

  variogramModel = list(model = model,params = pars)
  errSum = 0
  vario = data.matrix(varioIn)
  inp = which(names(varioIn) == "np")
  iacl1 = which(names(varioIn) == "acl1")
  iacl2 = which(names(varioIn) == "acl2")
  idist = which(names(varioIn) == "dist")
  igamma = which(names(varioIn) == "gamma")
  ia1 = which(names(varioIn) == "a1")
  ia2 = which(names(varioIn) == "a2")
  if (gDistEst) {
#dists is here the n*n matrix gDistObs
    if (dim(dists)[1] == dim(dists)[2]) {
      gd1 = mapply(FUN = function(i,dists) dists[i,i],vario[,iacl1],MoreArgs = list(dists))
      gd2 = mapply(FUN = function(i,dists) dists[i,i],vario[,iacl2],MoreArgs = list(dists))
      gb = mapply(FUN = function(i,j,dists) dists[i,j],vario[,iacl1],vario[,iacl2],MoreArgs = list(dists))
      gamma1 = mapply(FUN = varioEx,gd1,MoreArgs=list(variogramModel = variogramModel))
      gamma2 = mapply(FUN = varioEx,gd2,MoreArgs=list(variogramModel = variogramModel))
      gammab = mapply(FUN = varioEx,gb,MoreArgs=list(variogramModel = variogramModel))
      gammar = gammab-.5*(gamma1+gamma2)
      if (!is.null(aOver)) {
        farea = vario[,ia1]
        sarea = vario[,ia2]
        carea = mapply(FUN = function(i,j,aOver) aOver[i,j],vario[,iacl1],vario[,iacl2],MoreArgs = list(aOver))
        nugg = mapply(FUN = nuggEx,(1/farea + 1/sarea -2*carea/(farea*sarea))/2,
                      MoreArgs = list(variogramModel = variogramModel))
        gammar = gammar+nugg
      }
#         rnugg = (amp/farea+amp/sarea-2.*amp*aov/(farea*sarea))/2.
#        amp*(1/farea+1/sarea-2*aov/(farea+sarea))/2
    } else {
# Binned variograms - using geostatistical distance
# dists is here the n*3 matrix with distances within and between hypothetical areas
      gamma1 = mapply(FUN = varioEx,dists[,1],MoreArgs=list(variogramModel = variogramModel))
      gamma2 = mapply(FUN = varioEx,dists[,2],MoreArgs=list(variogramModel = variogramModel))
      gammab = mapply(FUN = varioEx,dists[,3],MoreArgs=list(variogramModel = variogramModel))
      gammar = gammab-.5*(gamma1+gamma2)
      if (!is.null(aOver)) {
        farea = vario[,ia1]
        sarea = vario[,ia2]
        carea = aOver
        nugg = mapply(FUN = nuggEx,(1/farea + 1/sarea -2*carea/(farea*sarea))/2,
                      MoreArgs = list(variogramModel = variogramModel))
        gammar = gammar+nugg
      }
    }
  } else {
# dists is here dObs
    if (inherits(varioIn,"rtopVariogramCloud")) {
      ar1 = mapply(FUN=function(i,dObs) dObs[[i]],vario[,iacl1],MoreArgs = list(dObs = dists))
      ar2 = mapply(FUN=function(i,dObs) dObs[[i]],vario[,iacl2],MoreArgs = list(dObs = dists))
      gammar = mapply(FUN = vred,a1 = ar1,a2 = ar2,MoreArgs = list(vredTyp = "ind",variogramModel = variogramModel))
      if (!is.null(aOver)) {
        farea = vario[,ia1]
        sarea = vario[,ia2]
        carea = mapply(FUN = function(i,j,aOver) aOver[i,j],vario[,iacl1],vario[,iacl2],MoreArgs = list(aOver))
        nugg = mapply(FUN = nuggEx,(1/farea + 1/sarea -2*carea/(farea*sarea))/2,
                      MoreArgs = list(variogramModel = variogramModel))
        gammar = gammar+nugg
      }
    } else {
# Binned variogram, not geostatistical distance
      gammar = mapply(FUN = vred,a1 = vario[,ia1],a2 = vario[,ia2],dist=vario[,idist],
                 MoreArgs = list(vredTyp = "hyp",variogramModel = variogramModel, resol = resol))
      if (!is.null(aOver)) {
        farea = vario[,ia1]
        sarea = vario[,ia2]
        carea = aOver
        nugg = mapply(FUN = nuggEx,(1/farea + 1/sarea -2*carea/(farea*sarea))/2,
                      MoreArgs = list(variogramModel = variogramModel))
        gammar = gammar+nugg
      }
    }
  }
  err = mapply(FUN=goFit,vario[,igamma],gammar,dist=vario[,idist],np = vario[,inp],MoreArgs=list(fit.method))
  errSum = sum(err)/sum(vario[,inp])
  if (debug.level > 0) print(paste(paste(round(pars,4), collapse = " "), "errSum = ", errSum))
  if (last) return(list(varFit = data.frame(vario,vario[,igamma], gammar = gammar, err = err),errSum = errSum)) else  return(errSum)
}





goFit = function(gobs,gest,dist,np,fit.method = 8) {
  if (fit.method == 1) {
    ww = np
  } else if (fit.method == 2) {
    ww = np*(gobs/gest-1)^2
    return(ww)
  } else if (fit.method == 6) {
    ww = 1
  } else if (fit.method == 7) {
    ww = 1/dist^2
  } else if (fit.method == 8) {
    ww = np*(gobs/gest-1)^2
    return(ww)
  } else if (fit.method == 9) {
    ww = np*(min((gobs/gest-1)^2,(gest/gobs-1)^2))
    return(ww)
  } else     stop(paste("fit.method", fit.method, "not implemented"))
  ww*(gobs-gest)^2
}


#cfunc = function(xd,param){
#  a = params[1]
#  b = params[2]
#  c = params[3]
#  d = params[4]

#  cfunc = a*(xd^b)*(1-exp(- ((xd/c)^d)))
#}





varioEx = function(skor,variogramModel) {
  model = variogramModel$model
  params = variogramModel$params
  res = 0.
  imod = imodel(model)
  vres = .Fortran("varioex",res,skor,length(params),params,imod)
#  print(paste(res,vres[[1]],params[1],params[2],params[3],params[4],params[5], sep = " " ))
  return(vres[[1]])
}

imodel = function(model) {
#     The numbers should match the numbers of the Fortran-function vario
  as.integer(switch(model, Exp = 1, Ex1 = 2, Gau = 3, Ga1 = 4, Gho = 5, Sph = 6, Sp1 = 7, Fra = 8) )
}

nuggEx = function(ared,variogramModel) {
  model = variogramModel$model
  params = variogramModel$params
  res = 0.
  
#  vres = .Fortran("nuggex",res,length(params),params,model)
  return(params[3]*ared)
}
