#.onLoad <- function(lib,pkg){
#    library.dynam("longitudinalData",pkg,lib)
#}


setGeneric("longData",function(traj,idAll,time,timeInData,varNames,maxNA){standardGeneric("longData")})
setGeneric("longData3d",function(traj,idAll,time,timeInData,varNames,maxNA){standardGeneric("longData3d")})
setGeneric("partition",function(clusters,traj,details=character()){standardGeneric("partition")})
setGeneric("imputation",function(traj,method="copyMean",lowerBound="globalMin",upperBound="globalMax"){standardGeneric("imputation")})
setGeneric("restoreRealData",function(object){standardGeneric("restoreRealData")})
setGeneric("expandParLongData",function(xParLongData,y){standardGeneric("expandParLongData")})
#setGeneric("plotTraj",function(x,y,...){standardGeneric("plotTraj")})
#setGeneric("plotSenators",function(x,y,...){standardGeneric("plotSenators")})
#setGeneric("plotMeans",function(x,y,...){standardGeneric("plotMeans")})
setGeneric("plotTraj3d",function(x,y,...){standardGeneric("plotTraj3d")})
setGeneric("plotMeans3d",function(x,y,...){standardGeneric("plotMeans3d")})
setGeneric("plot3dPdf",function(x,y,varY=1,varZ=2){standardGeneric("plot3dPdf")})
setGeneric("qualityCriterion",function(traj,clusters,imputationMethod="copyMean"){standardGeneric("qualityCriterion")})
#setGeneric("initializePartition",function(nbClusters,lengthPart,method="kmeans++",data){standardGeneric("initializePartition")})
setGeneric("plotCriterion",function(x,criterion=x["criterionActif"],nbCriterion=100){standardGeneric("plotCriterion")})
setGeneric("plotAllCriterion",function(x,criterion=CRITERION_NAMES[1:5],standardized = TRUE){standardGeneric("plotAllCriterion")})

#setGeneric("getClusters",function(cld, nbCluster, clusterRank = 1, asInteger = FALSE){standardGeneric("getClusters")})

#setGeneric("plotSubGroups",function(x,y,...){standardGeneric("plotSubGroups")})
#setGeneric("criterion",function(object,partition,method="linearInterpolation"){standardGeneric("criterion")})
#setGeneric("criterion3",function(object,partition,method="linearInterpolation"){standardGeneric("criterion3")})
#setGeneric("as.longData",function(data,id="I1",timeCol=0,timeReal=0,varName="V",other=list()){standardGeneric("as.longData")})

#setGeneric("selectSupTrajMinSize",function(object,minSize){standardGeneric("selectSupTrajMinSize")})
#setGeneric("expandPartition",function(object,listId){standardGeneric("expandPartition")})
#setGeneric("imputation",function(object,method,partition)standardGeneric("imputation"))
#setGeneric("plotTraj",function(x,y,...){standardGeneric("plotTraj")})
#setGeneric("plotAll",function(x,y,...){standardGeneric("plotAll")})




#setGenericVerif <- function(x,y){if(!isGeneric(x)){setGeneric(x,y)}else{}}
