#trajLong <- wideToLong(data.frame(id=as.character(1:200),traj),varying=2:11,id="id")
#names(trajLong)
#id="id"
#varying="t0"
#times="time"





plotShape <- function(toPlot,trajLong,centers,col="grey",colMeans,pourcent=0){
    if(toPlot=="traj"){
        matplotLong(trajLong,col=col,lwd=1,xlab="Times",main="",pourcent=pourcent)
    }else{
        if(toPlot=="means"){
            matplotLong(centers,col=1,lwd=8,lty=1,main="",pourcent=pourcent)
            matplotLong(centers,col=colMeans,lwd=4,lty=1,add=TRUE)
        }else{
            if(toPlot=="both"){
                matplotLong(trajLong,col=col,lwd=1,xlab="Times",main="",pourcent=pourcent)
                matplotLong(centers,col=1,lwd=8,add=TRUE,lty=1)
                matplotLong(centers,col=colMeans,lwd=4,add=TRUE,lty=1)
            }else{}
        }
    }
    return(invisible())
}



kmlShapeAux <- function (trajLong, seeds, FrechetSumOrMax, timeScale, aggregationMethod, shuffle, sampleSize, methodHclust, toPlot, maxIter){

    nbClusters <- length(seeds)
    id <- unique(trajLong[, 1])
    weight <- rep(1,length(id))
    centers <- trajLong[trajLong[, 1] %in% seeds,]
    centers[, 1] <- as.integer(factor(centers[, 1]))

    plotShape(toPlot=toPlot,trajLong=trajLong,centers=centers,colMeans=2:(nbClusters+1))

    iter <- 0
    group <- 0
    exGroup <- 1
    while (!identical(group, exGroup) && iter<maxIter ) {
        iter <- iter + 1
        exGroup <- group
        matriceDist <- numeric()
        for (iCenters in 1:nbClusters) {
            fDist <- function(i) {
                distFrechet(Px = trajLong[trajLong[, 1] == i,2],
                   Py = trajLong[trajLong[,1] == i, 3],
                   Qx = centers[centers[,1] == iCenters, 2],
                   Qy = centers[centers[,1] == iCenters, 3],
                   timeScale = timeScale,
                   FrechetSumOrMax = FrechetSumOrMax)
            }
            distToMi <- lapply(id, fDist)
            matriceDist <- cbind(matriceDist, distToMi)
        }
        group <- as.integer(as.factor(apply(matriceDist, 1, which.min)))
        names(group) <- id
        groupLong <- group[match(trajLong[,1],names(group))]

        nbClusters <- max(group)
        centers <- do.call(rbind,lapply(1:nbClusters,
            f<-function(iCenters) {
                 data.frame(iCenters,
                            meanFrechet(trajLong[groupLong==iCenters, ], timeScale = timeScale, aggregationMethod = aggregationMethod,
						    sampleSize=sampleSize,shuffle=shuffle, FrechetSumOrMax=FrechetSumOrMax,methodHclust = methodHclust)
		     )
		}
        ))

        pourcent <- sapply(1:nbClusters,function(x)sum(weight[group==x]))/sum(weight)
        plotShape(toPlot=toPlot,trajLong=trajLong,centers=centers,pourcent=pourcent,col=group+1,colMeans=2:(nbClusters+1))
    }
    if(iter==maxIter){warning("[kmlShape:kmlShape]: the maximum number of iteration has been reach, the algorithm did not converge\n")}else{}

    reOrder <- rank(-table(group),ties.method="first")
    group <- reOrder[group]

    centers$iCenters <- reOrder[centers$iCenters]
    centers <- centers[order(centers$iCenters,centers$times),]

    pourcent <- sapply(1:nbClusters,function(x)sum(weight[group==x]))/sum(weight)
    plotShape(toPlot=toPlot,trajLong=trajLong,centers=centers,pourcent=pourcent,col=group+1,colMeans=2:(nbClusters+1))

    names(group) <- id
    return(list(clusters = group, centers = centers))
}


#set.seed(6)
#myClds <- clds8
#nbClusters <- 2
#timeScale=0.1
#aggregationMethod="hierarchical"
#shuffle = TRUE
#sampleSize=NA
#methodHclust="average"
##toPlot = "both"
#maxIter=100
#FrechetSumOrMax = "max"




kmlShape <- function (myClds, nbClusters = 3, timeScale = 0.1, FrechetSumOrMax = "max", toPlot = "both", parAlgo=parKmlShape()){

    nameObject <- deparse(substitute(myClds))
    if(myClds["senatorsAvailable"]){
        trajLong <- myClds["senators"]
    }else{
        if(myClds["longAvailable"]){
            trajLong <- myClds["trajLong"]
        }else{
            trajLong <- reshapeWideToLong(cbind(myClds["id"],myClds["trajWide"]))
        }
    }

    id <- unique(trajLong[, 1])
    if(length(nbClusters)==1){seeds <- sample(id, nbClusters)}else{}

#    if(toPlot!="none"){matplotLong(trajLong,col="grey",lty=1)}else{}

    result <- kmlShapeAux(trajLong=trajLong, seeds=seeds, FrechetSumOrMax=FrechetSumOrMax, timeScale=timeScale, toPlot=toPlot,
                          aggregationMethod=parAlgo["aggregationMethod"], shuffle=parAlgo["shuffle"], sampleSize=parAlgo["sampleSize"],
                          methodHclust=parAlgo["methodHclust"],maxIter=parAlgo["maxIter"])


    myClds["kmlShape"] <- TRUE
    myClds["trajMeans"] <- result$centers
    if(myClds["senatorsAvailable"]){
        myClds["clustersSenators"] <- as.factor(result$clusters)
        names(myClds["clustersSenators"]) <- id
        myClds["clustersSenators"][myClds["mySenator"]$senators]
        myClds["clusters"] <- myClds["clustersSenators"][match(myClds["mySenator"]$senators, id)]
    }else{
        myClds["clusters"] <- as.factor(result$clusters)
    }
    names(myClds["clusters"]) <- myClds["id"]
    assign(nameObject, myClds, envir = parent.frame())
    return(myClds)
}


#plotTraj
#plotSenators
#plotMeans
#getMeans
#getClusters
#getSenators
#convertCldShapeToClassic
#convertCldClasicToShape
