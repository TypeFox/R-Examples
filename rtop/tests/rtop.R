set.seed(1501)
#-----------------------------
library(rtop)
library(rgdal)
options(error = recover)
  # Read directly from shape-files in data directory
rpath = system.file("extdata",package="rtop")
observations = readOGR(rpath, "observations")
predictionLocations = readOGR(rpath, "predictionLocations")
  #Finding a few prediction locations of them
  
  observations = observations[1:30,]
  predictionLocations = predictionLocations[1:2,]
  
  observations$obs = observations$QSUMMER_OB/observations$AREASQKM
  
  # Setting some parameters 
  params = list(gDist = TRUE, cloud = FALSE, rresol = 25, hresol = 3)
  # Build an object
  rtopObj = createRtopObject(observations,predictionLocations, params = params)
  # Fit a variogram (function also creates it)
  rtopObj = rtopFitVariogram(rtopObj)
  #rtopObj = checkVario(rtopObj)
  rtopObj$variogramModel                                                                        
  rtopObj2 = rtopKrige(rtopObj, cv = TRUE)


  print(attr(rtopObj2$varMatObs,"variogramModel"))
  
  rtopObj3 = rtopKrige(rtopObj)


 varmat = varMat(observations, predictionLocations, variogramModel = rtopObj$variogramModel, 
                 gDistEst = TRUE, gDistPred = TRUE, rresol = 25, hresol = 3)

all.equal(varmat$varMatObs, rtopObj2$varMatObs)
rtopObj4 = rtopKrige(rtopObj2)

#debug(rtop:::rtopDisc.SpatialPolygons)
  rtopObj5 = rtopKrige(rtopObj, params = list(cnAreas = 5, cDlim = 10, nclus = 2))
  
  print(summary(rtopObj2$predictions))
  print(summary(rtopObj3$predictions))
  print(summary(rtopObj4$predictions))
  print(all.equal(rtopObj4$predictions, rtopObj3$predictions))
  #spplot(rtopObj$predictions,col.regions = bpy.colors(), c("var1.pred","var1.var"))
  
  # Cross-validation
  #spplot(rtopObj2$predictions,col.regions = bpy.colors(), c("observed","var1.pred"))
  print(cor(rtopObj2$predictions$observed,rtopObj2$predictions$var1.pred))
  



  set.seed(1501)
  library(intamap)
  useRtopWithIntamap()
  library(intamap)
  output = interpolate(observations,predictionLocations,
     optList = list(formulaString = obs~1, gDist = TRUE, cloud = FALSE, nmax = 10, rresol = 25, hresol = 3), 
        methodName = "rtop")
  

print(all.equal(rtopObj4$predictions@data$var1.pred, output$predictions@data$var1.pred))
  print(all.equal(rtopObj4$predictions@data$var1.var, output$predictions@data$var1.var))


# Updating variogramModel
  
  rtopObj5 = varMat(rtopObj4)
  rtopObj6 = updateRtopVariogram(rtopObj5, exp = 1.5, action = "mult")
  rtopObj7 = varMat(rtopObj6)


rtopObj10 = rtopSim(rtopObj, nsim = 5)
rtopObj11 = rtopObj
rtopObj11$predictionLocations = rtopObj11$observations
rtopObj11$observations = NULL
rtopObj12 = rtopSim(rtopObj11, nsim = 10, beta = 0.01)

rtopObj10$simulations@data
rtopObj12$simulations@data


