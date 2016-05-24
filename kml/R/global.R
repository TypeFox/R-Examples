.onLoad <- function(lib,pkg){library.dynam("kml",pkg,lib)}

#setGeneric("clusterLongData",function(traj,idAll,time,timeInData,varNames,maxNA){standardGeneric("clusterLongData")})
#setGeneric("plotCriterion",function(x,criterion=CRITERION_NAMES[c(1,4,5)], standardized = TRUE){standardGeneric("plotCriterion")})
setGeneric("plotTraj",function(x,y,...){standardGeneric("plotTraj")})
setGeneric("plotMeans",function(x,y,...){standardGeneric("plotMeans")})
#setGeneric("kml",function(object,nbClusters=2:6,nbRedrawing=3,toPlot="none",parAlgo=parKml(),criterionNames=c("Calinski.Harabatz","test")){standardGeneric("kml")})
setGeneric("exportPartition",function(object,nbClusters,rank,nameObject,typeGraph="bmp",parTraj=parTRAJ(),parMean=parMEAN(),parWindows=windowsCut(1))standardGeneric("exportPartition"))
setGeneric("choice",function(object,typeGraph="bmp"){standardGeneric("choice")})
#setGeneric("plotA",function(x,y,parTraj=parTRAJ(),parMean=parMEAN(),parWin=windowsCut(x['nbVar']),nbSample=200,toPlot=c("both"),
 #                    criterion=x["criterionActif"],nbCriterion=100,standardized = FALSE)){standardGeneric("plotA")})
