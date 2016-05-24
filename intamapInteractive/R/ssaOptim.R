`ssaOptim` <-
function (observations, predGrid, candidates, action, nDiff,
    model, nr_iterations, plotOptim = TRUE, formulaString = NULL, 
    coolingFactor = nr_iterations/10, covariates = "over", ...)
{

  if (missing(formulaString) || is.null(formulaString)) {
    observations = SpatialPointsDataFrame(observations,
          data = data.frame(dum = rep(1, dim(coordinates(observations))[1])))
    formulaString = dum ~ 1
  }
    nn = dim(coordinates(observations))[1]
    nGrd = dim(coordinates(predGrid))[1]
    cnames = dimnames(coordinates(observations))[[2]]
    netPtsInit = observations
    lmcfit = NULL
    if (!is.null(formulaString) && length(attr(terms(formulaString), "term.labels")) > 0) {
      predGrid = SpatialPointsDataFrame(coordinates(predGrid),
                       data = model.frame(delete.response(terms(formulaString)),
                       predGrid), #proj4string = CRS(as.character(projs)), 
                       coords.nrs = c(1, 2))    
      observations = SpatialPointsDataFrame(coordinates(observations),
                       data = model.frame(terms(formulaString), observations), 
                       #proj4string = CRS(as.character(projs)),
                       coords.nrs = c(1, 2))
      dnp = which(names(predGrid) %in% dimnames(coordinates(predGrid))[[2]])
      if (length(dnp) > 0) { 
        predGrid = predGrid[, -dnp] #removes the matching xy predictors in predGrid
        if (dim(predGrid)[2] ==0) predGrid = SpatialPoints(predGrid)
      }
      dno = which(names(observations) %in% dimnames(coordinates(observations))[[2]])
      if (length(dno) > 0) { 
        observations = observations[, -dno] #removes the matching xy predictors in predGrid
        if (dim(observations)[2] ==0) observations = SpatialPoints(observations)
      }
      if (is.null(dim(predGrid))) covariates = NULL
    } else {
      predGrid = SpatialPoints(predGrid)
      observations = observations[,as.character(formulaString[[2]])]
      covariates = NULL
    }
    
   
    gridded(predGrid) = TRUE
    if (action == "del") { 
        sn = sample(nn)
        delPts = observations[which((sn) < (min(sn) + nDiff)),]
        netPts = observations[which((sn) >= (min(sn) + nDiff)),]
        addPts = NULL
        nn = dim(coordinates(netPts))[1]
        crit1 = calculateMukv(observations = netPts, predGrid = predGrid,
            model = model, formulaString = formulaString, ...)
        res = ssaMap(candidates, predGrid, model, max_points_shift = 1,
            maxShiftFactorX = 0.2, minShiftFactorX = 0, maxShiftFactorY = 0.2,
            minShiftFactorY = 0, start_p = 0.2, 
            netPts, addPts, delPts, crit1, nn, action, nDiff,
            netPtsInit, nr_iterations, plotOptim, formulaString = formulaString,
            coolingFactor = coolingFactor,
            ...)

        if ("data.frame" %in% getSlots(class(netPts)) && !"data.frame" %in%
            getSlots(class(res)))
            res = SpatialPointsDataFrame(res, 
                   data = netPts@data[as.numeric(row.names(res@data)), ])
        return(res)
    }
    if (action == "add") {
        cnames = dimnames(coordinates(observations))[[2]]
        delPts = NULL  # Just for the variable to exist when calling ssaMap
############ 
#        Added points are sampled from candidates
#        addPts = predGrid[which(sn <= nDiff), ]
#        netPts = rbind(observations, addPts)
############ added a proj4string to ensure matching projections when overlay is used
        projs = proj4string(predGrid)
        proj4string(candidates) = projs        
############ sample added pts from candidates
        addPts = spsample(candidates, nDiff, "random")

        if (!is.null(covariates)) {
          if (covariates == "krige") {
            dots = list(...)
            if ("nmax" %in% names(dots)) {
              nmax = dots$nmax
            } else nmax = 200
            addterms = as.data.frame(matrix(NA, ncol = length(names(observations)), 
                 nrow = nDiff))
            models = list()
            for (i in 1:length(names(observations))) {
              vname = names(observations)[i]
              names(addterms)[i] = vname
              if (names(observations)[i] == formulaString[[2]]) {
                models[[i]] = NULL
                addterms[,i] = 1 # Dummy for dependent variable
              } else {
                models[[i]] = autofitVariogram(as.formula(paste(names(observations)[i],
                        "~1")), predGrid)$var_model
                if (models[[i]]$model[2] == "Gau" & models[[i]]$psill[1] == 0)
                    models[[i]]$psill[1] = models[[i]]$psill[2]/100

                lres = krige(as.formula(paste(vname,"~1")), predGrid, 
                       addPts, model = models[[i]], nmax = nmax, debug.level = 0)$var1.pred
                if (is.factor(predGrid@data[,vname])) lres = factor(round(lres),
                   levels = levels(predGrid@data[,vname]))
                addterms[,i] = lres
              }
            }
          }  else {
            addterms = over(SpatialPoints(addPts, proj4string = CRS(as.character(projs))), predGrid)
### In case addpts have been sampled outside the area of observations
            gerr = predGrid@data[1,]
            for (i in 1:nDiff) if (is.na(addterms[i,1])) addterms[i,] = gerr
            addterms = cbind(addterms, rep(1,nDiff))
            names(addterms)[dim(addterms)[2]] = as.character(formulaString[[2]]) 
          }
############ addPts needs to be an SPDF to rbind with netPts
        } else addterms = 
            setNames(data.frame(rep(0, length(addPts))), as.character(formulaString[[2]]))
        if (sum(cnames %in% attr(terms(formulaString), "term.labels")) > 0)
            addterms = cbind(addterms, coordinates(addPts))
        addPts = SpatialPointsDataFrame(coordinates(addPts),
                       data = model.frame(terms(formulaString), addterms), 
                       proj4string = CRS(as.character(projs)),
                       coords.nrs = c(1, 2))
        dna = which(names(addPts) %in% dimnames(coordinates(observations))[[2]])
        if (length(dna) > 0) addPts = addPts[, -dna] #removes the matching xy predictors in predGrid
          
        netPts = rbind(observations, addPts)
        nn = dim(coordinates(netPts))[1]
        crit1 = calculateMukv(observations = netPts, predGrid = predGrid,
            model = model, formulaString = formulaString, ...)
        res = ssaMap(candidates, predGrid, model, max_points_shift = 1,
            maxShiftFactorX = 0.2, minShiftFactorX = 0, maxShiftFactorY = 0.2,
            minShiftFactorY = 0, start_p = 0.2, 
            netPts, addPts, delPts, crit1, nn, action, nDiff,
            netPtsInit, nr_iterations, plotOptim, formulaString = formulaString,
            models = models, coolingFactor = coolingFactor, covariates = covariates, ...)
        if ("data" %in% names(getSlots(class(netPts))) && 
            !"data" %in% names(getSlots(class(res))))
            res = SpatialPointsDataFrame(res, data = netPts@data[as.numeric(row.names(res@data)),
                ])
        return(res)
    }
}
