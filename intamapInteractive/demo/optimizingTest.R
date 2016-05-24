options(error = recover)
test = TRUE
set.seed(1)
library(intamapInteractive)
require(maptools)
# for SIC2004 dataset
data(sic2004)
coordinates(sic.val) = ~x+y
observations = sic.val["dayx"] 
coordinates(sic.grid)=~x+y
predGrid = sic.grid

#Finding the polygon for the candidate locations
bb = bbox(predGrid)
boun = SpatialPoints(data.frame(x=c(bb[1,1],bb[1,2],bb[1,2],bb[1,1],bb[1,1]),
                                y=c(bb[2,1],bb[2,1],bb[2,2],bb[2,2],bb[2,1])))
Srl = Polygons(list(Polygon(boun)),ID = as.character(1))
candidates = SpatialPolygonsDataFrame(SpatialPolygons(list(Srl)),
                                      data = data.frame(ID=1))

# Limits the number of prediction locations to have faster UK 
# computations
nGrid = dim(coordinates(predGrid))[1]
predGrid = predGrid[sample(seq(1,nGrid),1000),]
# Fits the variogram model (using function fit.variogram from package
# gstat)
model = fit.variogram(variogram(dayx~x+y, sic.val), vgm(50, "Sph", 250000, 250))
#plot(variogram(dayx~x+y, sic.val), model=model)
# Computes the Mukv of the current network
initMukv <- calculateMukv(observations, predGrid, model, formulaString = dayx~x+y)
print(initMukv)
# Computes Kriging and plot result 
GammaDoseMap = krige(dayx~x+y, observations, predGrid, model)
GammaDoseMap = as.data.frame(GammaDoseMap)
#windows()
levelplot(var1.pred~x+y, GammaDoseMap, aspect = "iso", col.regions=bpy.colors,
               panel = function(...) {
                       panel.levelplot(...)
                       panel.xyplot(y=sic.val$y, x=sic.val$x, col="white", pch=19);  panel.xyplot(y=sic.val$y, x=sic.val$x, col="black")
               }, main = "universal kriging prediction")
windows()
levelplot(var1.var~x+y, GammaDoseMap, aspect = "iso",
col.regions=bpy.colors,
               panel = function(...) {
                       panel.levelplot(...)
                       panel.xyplot(y=sic.val$y, x=sic.val$x, col="white", pch=19);  panel.xyplot(y=sic.val$y, x=sic.val$x, col="black")
                       
               }, main = "universal kriging variance")






###############################################################
# Deleting
###############################################################

# Deletes manually 20 stations from current network with method
# "manual" 
optimDel1=optimizeNetwork( observations,
                           method = "manual",
                           action = "del",
                           nDiff = 2,
                           predGrid, candidates, plotOptim = FALSE,  formulaString = dayx~x+y)
# Computes the Mukv of the optimized network with spatial # coverage
MukvDel1 <- calculateMukv(optimDel1, predGrid, model, formulaString = dayx~x+y)
print(MukvDel1)

# Deletes optimally 20 stations from current network with method
# "spcov" (spatial coverage)
optim1 = optimizeNetwork(observations,
         	                method = "spcov",
                          action = "del",
                          nDiff = 2,
                          predGrid, candidates,
                          plotOptim=TRUE)
# Computes the Mukv of the optimized network with spatial coverage
MukvDel2 = calculateMukv(optim1, predGrid, model, formulaString = dayx~x+y)
print(MukvDel2)

# Deletes optimally 20 stations from current network with method "ssa"
# (spatial simulated annealing) and criterion "mukv"
#windows()
optim2 = optimizeNetwork(observations ,
                          method = "ssa",
                          criterion = "MUKV",
                          action = "del",
                          nDiff = 2,
                          predGrid, candidates, model,
                          plotOptim=TRUE)
# Computes the Mukv of the optimized network with spatial simulated
# annealing applied to mukv
MukvDel3 <- calculateMukv(optim2, predGrid, model)
print(MukvDel3)

###############################################################
# Adding
###############################################################

# Adds manually 20 stations from current network with method
# "manual" 

optimAdd1=optimizeNetwork( observations,
                           method = "manual",
                           action = "add",
                           nDiff = 2,
                           predGrid, candidates)
# Computes the Mukv of the optimized network with spatial # coverage
MukvAdd1 <- calculateMukv(optimAdd1, predGrid, model)
print(MukvAdd1)


# Adds optimally 20 stations from current network with
# method "spcov" (spatial coverage)

  optimAdd2=optimizeNetwork( observations,
                           method = "spcov",
                           action = "add",
                           nDiff = 2,
                           predGrid, candidates, nGridCells = 5000,
                           nTry = 100, plotOptim=FALSE)
# Computes the Mukv of the optimized network with spatial # coverage
MukvAdd2 <- calculateMukv(optimAdd2, predGrid, model)
print(MukvAdd2)


# Adds optimally 20 stations from current network with
# method "ssa" (spatial simulated annealing) and
# criterion "mukv"

optimAdd3=optimizeNetwork( observations ,
                           method = "ssa",
                           criterion = "MUKV",
                           action = "add",
                           nDiff = 2,
                           predGrid, candidates, model,
                           plotOptim=FALSE, nr_iterations = )

# Computes the Mukv of the optimized network with spatial
# simulated annealing applied to mukv
MukvAdd3 <- calculateMukv(optimAdd3, predGrid, model)
print(MukvAdd3)


# Compares computed designs based on the Mukv results (Small MUKV is better)
print(initMukv)
# Deleting 20 measurements
print(MukvDel1) # Manual
print(MukvDel2) # Spatial Coverage
print(MukvDel3) # Spatial simulated annealing (MUKV in objective function)
# Adding 20 measurements
print(MukvAdd1) # Manual
print(MukvAdd2) # Spatial Coverage
print(MukvAdd3) # Spatial simulated annealing (MUKV in objective function)




