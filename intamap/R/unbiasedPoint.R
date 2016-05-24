

unbiasedKrige = function(object,formulaString,observations,predictionLocations,model,
     outputWhat,nmax,nsim,yamamoto,iwqmaxit = 500, iwqCpAddLim = 0.0001, debug.level,...) {
  dots = list(...)
  if (is(object,"Spatial")) {
    predictions = object
    iwqs = outputWhat[names(outputWhat) == "IWQSEL"]
    if (missing(debug.level)) debug.level = 1
    if (length(iwqs) > 0) {
      if (missing(formulaString)) {
        formulaString = as.formula(paste(names(observations)[1],"~1"))
        print(paste("warning: formulaString not given, using ",formulaString))
      }
      if (missing(nsim) || nsim == 0) nsim = 100
      if (missing(nmax)) nmax = 10  
    }
    if (missing(predictionLocations)) predictionLocations = predictions
  } else {
    outputWhat = object$outputWhat
    iwqs = outputWhat[names(outputWhat) == "IWQSEL"]
    predictions = object$predictions
    params = object$params
    if (missing(debug.level)) debug.level = object$debug.level
    if (length(iwqs) > 0) {
      model = object$variogramModel
      formulaString = object$formulaString
      observations = object$observations
      if (missing(nsim)) nsim = object$params$nsim
      if (missing(nmax)) nmax = object$params$nmax
      if (nsim == 0) nsim = 100
      if (missing(predictionLocations)) predictionLocations = object$predictionLocations
    }
  }

# Creating the accumulative distribution function
  acdf = acdfDef(predictions,...)
#  
# IWQSEL method requested
  if (length(iwqs) > 0) {
    if ("var1.pred" %in% names(predictions)) {
#    Simulations necessary
      if (inherits(object,"yamamoto") | ("yamamoto" %in% names(dots) && dots$yamamoto)) { 
        zPred = yamamotoKrige (formulaString,observations,predictionLocations,
          nsim=nsim,nmax = nmax,model = model,...) 
      } else {
        zPred = krige(formulaString,observations,predictionLocations,
                 nsim=nsim,model = model,nmax=nmax,...)
        print("Finished simulations")
      }
    } else {
      zPred = predictions
    } 
    for (iwq in 1:length(iwqs)) {
      iwqThresh = iwqs[[iwq]]
      pThresh = acdfFind(acdf,iwqThresh,inv=TRUE)
      iname = paste("IWQSEL",iwqThresh,sep="")
      iwqselRes = iwqsel(zPred@data,acdf,iwqThresh,pThresh,maxit = iwqmaxit, 
          cpAddLim = iwqCpAddLim, debug.level = debug.level, ...)
      predictions[[iname]] = iwqselRes$zEst 
    }
  }
  
  moks = outputWhat[names(outputWhat) == "MOK"]
  if (length(moks) > 0) {
#  Modified ordinary kriging prediction
    for (imok in 1:length(moks)) {
      mokThresh = moks[[imok]]
      pThresh = acdfFind(acdf,mokThresh,inv=TRUE)
      zp = quantile(predictions$var1.pred,pThresh)
       mname = paste("MOK",mokThresh,sep="")
      predictions[[mname]] = predictions$var1.pred + mokThresh-zp
    }
  }
      
  if (is(object,"Spatial")) predictions else {
    if ("predPoint" %in% names(object)) {
      object$predPoint = predictions
    } else {
      object$predictions = predictions
    } 
    return(object)
  }
}


acdfFindArr = function(FB,zArr,zMin,zInc) {
  ord = order(zArr)
  acdfArr = array(c(1:length(zArr)))
  nFB = dim(FB)[1]
  for (iz in 1:length(zArr)) {
    zm = zArr[ord[iz]]
    id =  as.integer(((zm-zMin)/zInc)+2)
    if (id >nFB) id = nFB
    if (id <1) id = 1
    acdfArr[ord[iz]] = FB[id]
  }     
  return(acdfArr) 
}



acdfFindSimp = function(FB,zVal,zMin,zInc) {
  if (missing(zMin) | missing(zInc)) {
    return(acdfFind(FB,zVal))
  } else {
#    iv =as.integer((zVal-zMin)/zInc)+2
#    tval = FB[iv,1]
#    p = FB[iv,2]
#    cat(paste(iv,zVal,tval,p,"\n"))
    id = as.integer(((zVal-zMin)/zInc)+2)
    if (id > dim(FB)[1]) return(1)
    if (id < 1) return(0)
    return(FB[id,2])
  }
}
# FBfindSimp(FBk,2.5,zMin = zmin,zInc = zinc)


acdfFind = function(FB,zVal,inv=TRUE,simp=FALSE,...) {
# Consider using findInterval
  if (simp & !inv) {
    if (missing(zMin) | missing(zInc)) {
      zMin = FB[[1,1]]
      zInc = FB[[2,1]]-zMin
    }
    id = as.integer(((zVal-zMin)/zInc)+2)
    if (id > dim(FB)[1]) return(1)
    if (id < 1) return(0)
    return(FB[id,2])
  } else {
    icol = ifelse(inv,1,2)
    jcol = ifelse(inv,2,1)
    il = max(which(FB[,icol] <= zVal))
    ih = min(which(FB[,icol] > zVal))
    b = FB[[ih,icol]]-FB[[il,icol]]
    b1 = (zVal-FB[[il,icol]])/b
    b2 = (FB[[ih,icol]]-zVal)/b
    if (il < 1) {
      il = 1
      b1 = 0
      b2 = 1
    }
    if (ih > dim(FB)[1]) {
      ih = dim(FB)[1]
      b1 = 1
      b2 = 0
    }
  }
#  cat(paste("FBfind",icol,jcol,zVal,il,ih,b,b1,b2,"\n"))
  return(b1*FB[ih,jcol] + b2*FB[il,jcol])
}



acdfDef = function(pPred,nfb = 200,...) {
# Creating the ACDF from ordinary kriging predictions
  FB = as.matrix(data.frame(ia = c(1:nfb),0))
  if ("var1.pred" %in% names(pPred)) {
  # kriging prediction
    nPred = dim(pPred@data)[1]
    zmin = min(pPred$var1.pred)-3*sqrt(max(pPred$var1.var))
    zmax = max(pPred$var1.pred)+3*sqrt(max(pPred$var1.var))
    zinc = (zmax-zmin)/(nfb-2)
    z0 = zmin-0.5*zinc
    for (ia in 1:nfb) {
      zVal = z0+zinc*(ia-1)
      FB[ia,1] = zVal
      FB[ia,2] = sum(pnorm((zVal-pPred$var1.pred)/sqrt(pPred$var1.var+0.00001)))/nPred
#    cat(paste(ia,round(zinc,3),round(zVal,3),round(FB[ia,1],3),round(FB[ia,2],6),"\n"))
    }
  } else if ("sim1" %in% names(pPred)) {
  # simulations
    zmin = min(pPred)
    zmax = max(pPred)
    zinc = (zmax-zmin)/(nfb-2)
    z0 = zmin-0.5*zinc
    pall = dim(pPred)[1] * dim(pPred)[2]
    for (ia in 1:nfb) {
      zVal = z0+zinc*(ia-1)
      FB[ia,1] = zVal
      FB[ia,2] = sum(pPred<zVal)/pall
    }
  }
  return(FB)
}




