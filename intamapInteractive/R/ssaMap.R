############################################################################
# Spatial simulated annealling with one criterion. 
# This function calls on 1 related functions (calculateMukv) 
# and returns criterionIterf (criterion for all iterations)
# Stopping criterion is given by a number of times 
# without improvement in search of a better design (countMax)
############################################################################

ssaMap = function(candidates, predGrid, model, max_points_shift, maxShiftFactorX, minShiftFactorX, 
                  maxShiftFactorY, minShiftFactorY, start_p,  
                  netPts, addPts, delPts, crit1, nn, action, nDiff, netPtsInit, 
                  nr_iterations, plotOptim,
                  formulaString = NULL, models, nmax = 200, coolingFactor = nr_iterations/10, 
                  covariates = "over", countMax = 200, ...) {
  
  
 
  cnames = dimnames(coordinates(netPts))[[2]]
# which are the coordinate variables of the data.frame of the locations  
  cvar = which(names(as.data.frame(netPts)) %in% cnames)

  nvar = dim(as.data.frame(netPts))[2]-2
  cform = as.formula(paste("~",cnames[1],"+",cnames[2]))
  count = 0
#### Using the values of first record in case points are sampled outside predGrid
  if (is.null(formulaString) || length(attr(terms(formulaString), "term.labels")) == 0 ||
     all(attr(terms(formulaString), "term.labels") %in% cnames)) 
    covariates = NULL
  if (!is.null(covariates) && covariates == "over") gerr = predGrid@data[1,]

# 
  if (missing(candidates)) {
    bb = bbox(netPts)
    bb1 = bbox(predGrid)
    bb[[1]] = min(bb[[1]],bb1[[1]])
    bb[[2]] = min(bb[[2]],bb1[[2]])
    bb[[3]] = max(bb[[3]],bb1[[3]])
    bb[[4]] = max(bb[[4]],bb1[[4]])
    boun = SpatialPoints(data.frame(x=c(bb[1,1],bb[1,2],bb[1,2],bb[1,1],bb[1,1]),
                                y=c(bb[2,1],bb[2,1],bb[2,2],bb[2,2],bb[2,1])))
    Srl = Polygons(list(Polygon(boun)),ID = as.character(1))
    candidates = SpatialPolygonsDataFrame(SpatialPolygons(list(Srl)),
                                      data = data.frame(ID=1))
  }
  
# settings for simulated annealing
  x_bounds <- bbox(candidates)[1, ]
  y_bounds <- bbox(candidates)[2, ]
  x_extent <- x_bounds[2] - x_bounds[1]
  y_extent <- y_bounds[2] - y_bounds[1]
  max_shift_x <- maxShiftFactorX * x_extent # maximum shift in x-direction <Jan-Willem hoger>
  min_shift_x <- minShiftFactorX * x_extent  # minimum shift in x-direction
  max_shift_y <- maxShiftFactorY * y_extent  # maximum shift in y-direction <Jan-Willem hoger, to 0.50>
  min_shift_y <- minShiftFactorY * y_extent  # minimum shift in y-direction
  
  nr_designs <- 1  # counter for number of accepted designs

  oldpoints <- netPts  # save initial design because new design may be rejected
  criterionInitial <- crit1
  oldcriterion <- criterionInitial  # also save current criterion
  oldDelPoints <- delPts
  criterionIterf <- NULL
  bestCriterion = Inf
  for (k in 1:nr_iterations) {
# scenario of deletion
    if (action == "del") {    

      selected_shifts_del <- sample(nDiff)  
      selected_shifts_net <- sample(nn)
      oldDelPt = oldDelPoints[which(selected_shifts_del<=max_points_shift),] 
      newDelPt = oldpoints[which(selected_shifts_net<=max_points_shift),]

      delPts=rbind(oldDelPoints[which(selected_shifts_del>max_points_shift),],newDelPt)
      netPts=rbind(oldpoints[which(selected_shifts_net>max_points_shift),],oldDelPt)
      criterion = calculateMukv(observations = netPts, predGrid = predGrid, model = model, 
                  formulaString = formulaString, ...)

    } else if (action=="add") {
# scenario of addition
      oldpoints = as.data.frame(oldpoints)
      netPts = as.data.frame(netPts)
      selected_shifts <- c(rep(max_points_shift+1,nn-nDiff),sample(seq(1,nDiff)))
      arraypos <- which(selected_shifts<=max_points_shift)
      selected_shifts[-arraypos] = 0
      outside = "TRUE"
      while (outside){
        x_shift <- max_shift_x-k/nr_iterations*(max_shift_x-min_shift_x)  
# possible shift in x-direction decreases linearly
        y_shift <- max_shift_y-k/nr_iterations*(max_shift_y-min_shift_y)  
# possible shift in y-direction decreases linearly
        netPts[,cvar[1]] <- oldpoints[,cvar[1]] + x_shift*
            ((as.matrix(rep(-1,nn))+2*as.matrix(runif(nn)))*as.matrix(selected_shifts))  
# true shift in x-direction
        netPts[,cvar[2]] <- oldpoints[,cvar[2]] + y_shift*
            ((as.matrix(rep(-1,nn))+2*as.matrix(runif(nn)))*as.matrix(selected_shifts))  
# true shift in y_direction   
        net=netPts
        coordinates(net)= cform
        while(length(zerodist(net))[1]>0){
          netPts[,cvar[1]] <- oldpoints[,cvar[1]] + x_shift*
              ((as.matrix(rep(-1,nn))+2*as.matrix(runif(nn)))*as.matrix(selected_shifts))  
# true shift in x-direction
          netPts[,cvar[2]] <- oldpoints[,cvar[2]] + y_shift*
              ((as.matrix(rep(-1,nn))+2*as.matrix(runif(nn)))*as.matrix(selected_shifts))  
# true shift in y_direction
          net=netPts
          coordinates(net)=cform
          }
        pointx <- netPts[arraypos, cvar[1]]
        pointy <- netPts[arraypos, cvar[2]]
        pointxy <- netPts[arraypos,]
        projxy = proj4string(candidates)
        outside <- is.na(over(SpatialPoints(matrix(c(pointx, pointy), ncol = 2), 
              proj4string = CRS(as.character(projxy))), candidates)[1])
      }
      
      coordinates(netPts) = cform
      if (!is.null(covariates)) {
        gridded(predGrid) = TRUE
        if (covariates == "over") {
           overselect<-over(SpatialPoints(netPts[arraypos,]), predGrid)
          if (is.na(overselect[1])) overselect = gerr
          for (i in 2:nvar) {
            if (is.factor(netPts@data[,i]) & !is.factor(overselect[i-1])) 
                overselect[i-1] = factor(overselect[i-1],levels = levels(netPts@data[,i]))
            netPts@data[arraypos,i] = overselect[i-1]
          }
        } else if (covariates == "krige") {
          netPts@data[arraypos,1] = 0
          for (i in 2:nvar) {
            lres = krige(as.formula(paste(names(netPts)[i],"~1")), predGrid, 
                          netPts[arraypos,], model = models[[i]], nmax = nmax, 
                          debug.level = 0)$var1.pred
            if (is.factor(netPts@data[,i])) lres = factor(round(lres),levels = levels(netPts@data[,i]))
            netPts@data[arraypos,i] = lres
          }
        } else stop(paste("Not able to use method", covariates, "for interpolating covariates"))
      }
      criterion = calculateMukv(observations = netPts, predGrid = predGrid, model = model, 
            formulaString = formulaString, ...)
      netPts <- as.data.frame(netPts) # need as dataframe for oldpoints
    }

    p = runif(1) # to allow accepting an inferior design
    if (criterion <= oldcriterion){
      oldpoints = netPts
      if (action == "del") oldDelPoints = delPts
      oldcriterion = criterion
      nr_designs = nr_designs+1
      count = 0
      cat("No improvement for",count,"iterations ")
      cat("[Will stop at",countMax,"iterations with no improvement]", "\n")
    } else if (criterion > oldcriterion & p <= (start_p*exp(-k/coolingFactor))){
      oldpoints = netPts
      if (action == "del") oldDelPoints = delPts
      oldcriterion = criterion
      nr_designs = nr_designs+1
      count = count + 1
      cat("No improvement for",count,"iterations  p = ", p, "lim = ",start_p*exp(-k/coolingFactor) )
      cat("[Will stop at",countMax,"iterations with no improvement]", "\n")
    } else{
      criterion = oldcriterion
      if (action == "del") oldDelPt = oldDelPoints
      netPts = oldpoints
      nr_designs = nr_designs
      count = count + 1
      cat("No improvement for",count,"iterations ")
      cat("[Will stop at",countMax,"iterations with no improvement]", "\n")
    }
    criterionIterf[k] <- criterion
    if (criterion < bestCriterion/1.0000001) {
      bestNetPts = netPts
      bestCriterion = criterion
      bestOldcriterion = oldcriterion
      bestOldpoints = oldpoints
      if (action == "del") {
        bestOldDelPoints = oldDelPoints
        bestOldDelPt = oldDelPt
        bestDelPts = delPts
      }
    }

    if (plotOptim){
      if (!missing(candidates)) {                     
        plot(candidates)
        if (action == "del") {
          points(oldpoints, col = 1, pch = 19, cex = 0.7)
          points(oldDelPoints, col = 2, pch = "X", cex = 1.2)
        } else {
          points(netPtsInit, col=1, pch = 19, cex = 0.7)  
          points(netPts[(nn-nDiff+1):nn,cvar],  col = "green", pch = 19)
          points(pointxy,  col = 2, pch = 19)
        }
      } else {
        plot(oldpoints, col = 1, pch = 19, cex = 0.7)
        points(oldDelPoints, col = 2, pch = "X", cex = 1.2)
      }
      title('Spatial Simulated Annealing', 
                xlab=(paste("Criterion = ",  signif(criterion, digits=5),
                       "(Best Criterion = ",  signif(bestCriterion, digits=5),")")), 
                ylab=(paste("Iterations = ", k)))
    }

    if (count == countMax) {
      if (criterion > bestCriterion*1.000001) {
        oldpoints = bestOldpoints
        netPts = bestNetPts
        oldcriterion = bestOldcriterion
        criterion = bestCriterion
        if (action == "del") {
          oldDelPoints = bestOldDelPoints 
          oldDelPt = bestOldDelPt 
          delPts = bestDelPts
        }
        nr_designs = nr_designs + 1
        count = 0
        cat("Reached countMax with suboptimal design, restarting with previously best design \n")
        cat("No improvement for",count,"iterations ")
        cat("[Will stop at",countMax,"iterations with no improvement]", "\n")
      } else break
    }
  } 

  if (!inherits(netPts,"Spatial")) coordinates(netPts) = cform
  attr(netPts,"criterion") = criterionIterf
  return(netPts)
}


  
  
