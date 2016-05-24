
rtopKrige.STSDF = function(object, predictionLocations = NULL, varMatObs, varMatPredObs,
                   varMat, params = list(), formulaString, sel,  
                   olags = NULL, plags = NULL, lagExact = TRUE, ...) {
  if (!requireNamespace("spacetime")) stop("spacetime not available")
  params = getRtopParams(params, ...)
  cv = params$cv
  debug.level = params$debug.level
  if (!cv & !all.equal(is.null(olags), is.null(plags))) 
    stop("Lag times have to be given for both observations and predictionLocations, or for none of them")
  if (!missing(varMat) && missing(varMatObs)) {
    if (is.atomic(varMat)) {
      varMatObs = varMat
      if (!cv) stop(paste("Not cross-validation, you must provide either varMatObs",
                          " and varMatPredObs or varMat with the matrices as elements"))
    } else {
      varMatObs = varMat$varMatObs
      varMatPredObs = varMat$varMatPredObs
    }
  }
  
  if (cv) {
    predictionLocations = object
    plags = olags

    varMatPredObs = varMatObs
  }
  ntime = dim(object)[2]
  nspace = dim(object)[1]
  indx = predictionLocations@index
  sinds = sort(unique(indx[,1]))
  tinds = sort(unique(indx[,2]))
  pspace = dim(predictionLocations)[1]
  ptime = dim(predictionLocations)[2]
  pfull = if (prod(dim(predictionLocations)) == dim(predictionLocations@data)[1]) TRUE else FALSE
  if (is(predictionLocations, "STSDF")) {
    predictionLocations@data = cbind(predictionLocations@data, 
                                   data.frame(var1.pred = NA, var1.var = NA, var1.yam = NA))
  } else if (is(predictionLocations, "STS")) {
    predictionLocations = spacetime::STSDF(predictionLocations, data = data.frame(var1.pred = NA, var1.var = NA, var1.yam = NA))
  }
  if (is.null(olags) & prod(dim(object)) == dim(object@data)[1]) { # if all stations have obs from all time steps
    obj1 = object[,1]
    if (is(predictionLocations, "STSDF")) {
      predLoc = predictionLocations@sp
    } else {
      predLoc = predictionLocations
    }
    ret = rtopKrige.default(obj1, predLoc, varMatObs,
                                   varMatPredObs, varMat, params, formulaString, wret = TRUE)#,
    #sel, ...)
    weight = ret$weight
    wvar = ret$predictions$var1.var
    obs = as.data.frame(object)
    obs = obs[,c("sp.ID", "timeIndex", as.character(formulaString[[2]]))]
 #   obs$timeIndex = object@index[,2]
    obs = reshape(obs, v.names = as.character(formulaString[[2]]),
                  idvar = "sp.ID", timevar = "timeIndex", direction = "wide")
    if (interactive() & debug.level == 1 & length(sinds) > 1) pb = txtProgressBar(1, length(sinds), style = 3)
    for (istat in 1:length(sinds)) {
      isind = sinds[istat]
      ttinds = which(indx[,1] == isind)
      lweight = weight[istat,]
      preds = lweight %*% as.matrix(obs[,2:dim(obs)[2]])
      diffs = sweep(obs[,2:dim(obs)[2]], 2, preds )
      var1.yam = t(lweight) %*% (diffs^2)
      #      var1.var = t(lweights) %*% ((as.matrix(Obs[iobs,depVar]@data)-as.numeric(var1.pred))^2)
      
      if (!pfull) {
        ptinds = indx[ttinds,1]
        preds = preds[ptinds]
        var1.yam = var1.yam[ptinds]
      }
      predictionLocations@data$var1.pred[ttinds] = preds
      predictionLocations@data$var1.var[ttinds] = wvar[istat]
      predictionLocations@data$var1.yam[ttinds] = var1.yam
      if (interactive() & debug.level == 1 & length(sinds) > 1) setTxtProgressBar(pb, istat)
    }
    if (interactive() & debug.level == 1 & length(sinds) >  1) close(pb)
  } else { 
    object@sp$sindex = sindex = 1:nspace
    object@time$tindex = tindex = 1:ntime
    predictionLocations@sp$sindex = 1:pspace
    predictionLocations@time$tindex = 1:ptime
    obs = as.data.frame(object)
    depVar = as.character(formulaString[[2]])
    obs = obs[,c("sp.ID", "timeIndex", depVar)]
    #   obs$timeIndex = object@index[,2]
    if (!requireNamespace("reshape2")) stop("reshape2 not available")
    obs = reshape2::dcast(obs, sp.ID ~ timeIndex)
    #        reshape(obs, v.names = as.character(formulaString[[2]]),
    #                    idvar = "sp.ID", timevar = "timeIndex", direction = "wide")
    obs = obs[order(as.numeric(as.character(obs$sp.ID))), ]
    obs = as.matrix(obs[,2:(ntime+1)])
    
    if (!cv | !is.null(olags)) {
      ploc = as.data.frame(predictionLocations)
      ploc = cbind(ploc, var = 1)
      ploc = ploc[,c("sp.ID", "timeIndex", "var")]
      ploc$sp.ID = as.numeric(as.character(ploc$sp.ID))
      ploc = reshape2::dcast(ploc, sp.ID ~ timeIndex)
      ploc = as.matrix(ploc)
    #    spobs = NULL
    }
    if (interactive() & debug.level == 1 & ntime >  1) pb = txtProgressBar(1, ntime, style = 3)
    spPredLoc = if (!cv) predictionLocations@sp else NULL
    if (is.null(olags)) {
      oldind = NULL
      newkrige = 0
      for (itime in 1:ptime) {
        newind = which(!is.na(obs[,itime]))
        if (length(newind) == 0) next
        #    print(paste(1, istat, itime, length(newind)))
        if (is.null(oldind) || !isTRUE(all.equal(newind, oldind))) {
          oldind = newind
          newkrige = newkrige + 1
          ppq = object[,itime]
          vmat = varMatObs[newind, newind]
          if (cv) {
            vpredobs = NULL 
          } else {       
            vpredobs = varMatPredObs[newind,]
          }
          ret = rtopKrige.default(object = ppq, predictionLocations = spPredLoc, varMatObs = vmat,
                                         varMatPredObs = vpredobs, params = params, 
                                         formulaString = formulaString, wret = TRUE, debug.level = 0)#,
          weight = ret$weight
          wvar = ret$predictions$var1.var
        }
        ob = obs[,itime]
        ob = ob[!is.na(ob)]
        preds = weight %*% ob
        
        var1.yam = rep(NA, dim(weight)[1])
        for (istat in 1:dim(weight)[1]) {
          diffs = ob - preds[istat]
          var1.yam[istat] = sum(weight[istat,] * t(diffs^2))  
        }
        
        itind = tinds[itime]
        ttinds = which(indx[,2] == itind)
        
        if (!pfull & FALSE) {
          ptinds = indx[ttinds,1]
          preds = preds[ptinds]
          var1.yam = var1.yam[ptinds]
        }
        
        predictionLocations@data$var1.pred[ttinds] = preds
        predictionLocations@data$var1.var[ttinds] = wvar
        predictionLocations@data$var1.yam[ttinds] = var1.yam
        if (interactive() & debug.level == 1 & ntime > 1) setTxtProgressBar(pb, itime)
      }
    } else {
      ntot = dim(predictionLocations@data)[1]
      if (interactive() & debug.level == 1 & ntime >  1) pb = txtProgressBar(1, ntot, style = 3)
      itot = 0
      objsp = object@sp
      objsp = SpatialPointsDataFrame(SpatialPoints(objsp), data = objsp@data)
      newkrige = 0
      obs2 = obs
      for (istat in 1:pspace) {
 #       neigh = which(distm < maxdist)
        
        ispace = predictionLocations@sp@data$sindex[istat]
        pxts = as.data.frame(predictionLocations[istat,])
        stpred = predictionLocations@sp[istat,]
        distm = spDistsN1(coordinates(object@sp),coordinates(stpred))
        ppred = SpatialPoints(stpred) 
        rolags = olags - plags[istat] # longer lags will be positive, i.e., use future observations
        nbefore = -floor(min(rolags))
        nafter = ceiling(max(c(rolags, 1)))
        naadd1 = matrix(NA, nrow = nspace, ncol = nbefore)#, dimnames = list(as.character(1:nbefore), names(obs)))
        naadd2 = matrix(NA, nrow = nspace, ncol = nafter)#, dimnames = list(names(obs)as.character(1:nafter)))
        obsb = cbind(naadd1, cbind(obs2, naadd2))
        obs[,] = NA
        if (!lagExact) {
          rolags = round(rolags, 0)
          for (jstat in 1:nspace) {
            obs[jstat, ] = obsb[jstat, (nbefore + 1 + rolags[jstat]) : (nbefore + ntime + rolags[jstat])]
          }
        } else {
          rdiff = rolags - floor(rolags)
          for (jstat in 1:nspace) {
            nf = (nbefore + 1 + floor(rolags[jstat])) : (nbefore + ntime + floor(rolags[jstat]))
            nl = (nbefore + 1 + ceiling(rolags[jstat])) : (nbefore + ntime + ceiling(rolags[jstat]))
            obs[jstat, ] = obsb[jstat, nf] * (1-rdiff[jstat]) + obsb[jstat, nl] * rdiff[jstat]
          }
        }
        oldind = NULL
        ntl = length(pxts$tindex)
        preds = rep(NA, ntl)
        var1.yam = rep(NA, ntl)
        var1.var = rep(NA, ntl)
        tms = as.numeric(as.character(pxts$tindex))
        if (cv) vorder = order(varMatObs[istat,]) else vorder = order(varMatPredObs[istat])
        for (itime in seq_along(tms)) {
          jtime = tms[itime]
          tp = ploc[istat, jtime]
          if (is.na(tp)) next
          newind = vorder[vorder %in% which(!is.na(obs[,jtime]))]
          if (cv) newind = newind[-which(newind == istat)]
          if (is.numeric(params$nmax) && length(newind) > params$nmax) newind = newind[1:params$nmax]
          if (length(newind) == 0) next
           
          newind = newind[newind %in% vorder]          
          #    print(paste(1, istat, itime, length(newind)))
          if (is.null(oldind) || !isTRUE(all.equal(newind, oldind))) {
            oldind = newind
            newkrige = newkrige + 1
      #      print(paste(2, istat, itime, length(newind), newkrige))
            ppq = objsp[newind, ]
            ppq$intvar = obs[newind, jtime]
            vmat = varMatObs[newind, newind]
            vpredobs = varMatPredObs[newind, istat, drop = FALSE]
        #    fmstring = as.formula(paste(names(ppq)[1], "~1"))
            ret = rtopKrige.default(object = ppq, ppred, varMatObs = vmat,
                                           varMatPredObs = vpredobs, params = params, 
                                           formulaString = intvar~1, wret = TRUE, 
                                           debug.level = 0, cv = FALSE)#,
            weight = ret$weight
            wvar = ret$predictions$var1.var
          }
          preds[itime] = sum(weight * obs[newind,jtime])
          diffs = obs[newind,jtime] - preds[itime]
          var1.yam[itime] = sum(weight * (diffs^2))
          var1.var[itime] = wvar
          itot = itot + 1
          if (interactive() & debug.level == 1 & ntime > 1) setTxtProgressBar(pb, itot)
        }
#        print(paste(1, istat, itime, length(newind), newkrige))
      
        itind = tinds[itime]
        ttinds = which(indx[,1] == ispace)
        
        predictionLocations@data$var1.pred[ttinds] = preds
        predictionLocations@data$var1.var[ttinds] = var1.var
        predictionLocations@data$var1.yam[ttinds] = var1.yam
      }
    }
    if (interactive() & debug.level == 1 & ntime >  1) close(pb)
  }
  list(predictions = predictionLocations)
}

