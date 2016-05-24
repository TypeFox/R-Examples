cat("\n####################################################################
############################### kml3d ##############################
############################# Creation #############################
####################################################################\n")


parKml3d <- function (saveFreq = 100, maxIt = 200, imputationMethod = "copyMean",
    distanceName = "euclidean3d", power = 2, distance = function() {
    }, centerMethod = meanNA, startingCond = "nearlyAll", nbCriterion = 100,scale=TRUE)
{
    if (distanceName == "euclidean3d") {
        distance <- function(x,y){dist(rbind(c(x),c(y)),method='euclidean')}
    }else{}
    new("ParKml", saveFreq = saveFreq, maxIt = maxIt, imputationMethod = imputationMethod,
        distanceName = distanceName, power = power, distance = distance,
        centerMethod = centerMethod, startingCond = startingCond,
        nbCriterion = nbCriterion, scale=scale)
}




### ATTENTION : calculCenterGeneralized travaillait avec traj + Partition ;
### Maintenant, il travaille avec traj + part.
#calculCenterGeneralized <- calculTrajMean <- function(traj,part,centerMethod=meanNA){
#    trajMean <- apply(traj, c(2,3), tapply, part, centerMethod)
#    return(trajMean)
#}

### ATTENTION : Ne fonctionne pas avec les partitions à un seul clusters
calculTrajMean3d <- function(traj,clust,centerMethod=function(x){mean(x,na.rm=TRUE)}){
    trajMean <- apply(traj, c(2,3), tapply, clust, centerMethod)
    return(trajMean)
}

### CalculMean : même chose mais en C.


### On suppose que si un centre est NA, il est en dernière ligne de clustersCenter
affectIndiv3d <- function(traj,clustersCenter,distance=dist3d){
#    if (distance %in% METHODS){distanceFun <- ,method=distance))}}else{distanceFun <- distance}
    nbId <- nrow(traj)
    clusterAffectation <- rep(1,nbId)
    distActuel <- apply(traj,1,function(x){distance(x,clustersCenter[1,,])})
    ##   print(distActuel)
    for(iNbClusters in 2:nrow(clustersCenter)){
        distToMean <- apply(traj,1,function(x){distance(x,clustersCenter[iNbClusters,,])})
 #       print(distToMean)
        cond <- distToMean<distActuel
        cond[is.na(cond)] <- FALSE # Car si cond==NA, c'est que distToMean==NA et donc on ne change pas l'affectation.
        clusterAffectation <- ifelse(cond,rep(iNbClusters,nbId),clusterAffectation)
        distActuel <- ifelse(distToMean<distActuel,distToMean,distActuel)
    }
    return(clusterAffectation)
}


kml3dSlow <- function(traj,clusterAffectation,toPlot="traj",parAlgo=parKml3d()){
#    if (distance %in% METHODS){distanceFun <- function(x,y){return(dist(t(cbind(x,y)),method=distance))}}else{distanceFun <- distance}
 #   print(distanceFun)
    longDat3dTraj <- longData3d(traj,maxNA=ncol(traj)-1)
    kmlCenterMethod <- parAlgo['centerMethod']
    kmlDistance <- parAlgo['distance']

    exClusterAffectation <- partition()
    if(toPlot%in%c("traj","both")){
        screens <- ClusterLongData3d_plotTrajMeans(longDat3dTraj,partition(clusterAffectation),addLegend=TRUE)
        close.screen(screens)
    }else{}
    for(iterations in 1:parAlgo['maxIt']){
        clustersCenter <- calculTrajMean3d(traj=traj,clust=clusterAffectation,centerMethod=kmlCenterMethod)
        clusterAffectation <- affectIndiv3d(traj=traj,clustersCenter=clustersCenter,distance=kmlDistance)
        if(toPlot%in%c("traj","both")){
           screens <- ClusterLongData3d_plotTrajMeans(longDat3dTraj,partition(clusterAffectation),addLegend=TRUE)
           close.screen(screens)
        }else{}
        if(identical(clusterAffectation,exClusterAffectation)){
            clusterAffectation <- partition(clusterAffectation,longDat3dTraj,
               details=c(convergenceTime=as.character(iterations),algorithm="kmeans 3d, slow (R)",multiplicity="1"))
            return(clusterAffectation)
        }else{
            exClusterAffectation <- clusterAffectation
        }
    }
    return(partition(clusterAffectation,longDat3dTraj,details=c(convergenceTime=as.character(Inf),algorithm="kmeans 3d, slow (R)",multiplicity="1")))
}


kml3dFast <- function(traj,clusterAffectation){
    traj <- matrix(traj,nrow=dim(traj)[1])
    ### Modifier le nom de l'algorithm qui trouve
    return(kmlFast(traj,clusterAffectation))
}


fastOrSlow3d <- function(toPlot,distName){
    if(toPlot%in%c("both","traj") | !(distName=="euclidean3d")){
        cat(" ~ Slow KmL3D ~\n")
        fast <- FALSE
    }else{
        cat(" ~ Fast KmL3D ~\n")
        fast <- TRUE
    }
    return(fast)
}


kml3d <- function(object,nbClusters=2:6,nbRedrawing=20,toPlot="none",parAlgo=parKml3d()){
    if(class(object)=="ClusterLongData"){
        stop("[kml3d]: kml3d is for joint longitudinal data (object 'ClusterLongData3d').
For classic longitudinal data (object of class 'ClusterLongData'), use kml")
    }else{}

    nameObject<-deparse(substitute(object))
    on.exit(if(toPlot!="none"){close.screen(listScreen)}else{})

    if(parAlgo["scale"]){scale(object)}else{}

    nbIdFewNA <- object["nbIdFewNA"]
    convergenceTime <- 0
    traj <- object["traj"]
    nbTime <- length(object["time"])
    saveCld <-0

    ################
    ## listScreen[1] (à droite) est pour les traj.
    listScreen <- cutScreen(toPlot)
    if(toPlot%in%c("both","criterion")){
        screen(listScreen[2])
        plotCriterion(as(object,"ListPartition"),nbCriterion=parAlgo['nbCriterion'])
    }else{}

    ################
    ## Starting conditions
    startingCond <- expandStartingCond(parAlgo['startingCond'],nbRedrawing,object["initializationMethod"])
    object["initializationMethod"] <- unique(c(object["initializationMethod"],startingCond))

    ################
    ## Fast or Slow, according to distance and to toPlot
    fast <- fastOrSlow3d(toPlot,parAlgo['distanceName'])

    for(iRedraw in 1:nbRedrawing){
        for(iNbClusters in nbClusters){
            saveCld <- saveCld+1
            clustersInit <- initializePartition(nbClusters=iNbClusters,lengthPart=nbIdFewNA,method=startingCond[iRedraw],data=traj)
            clust <- rep(NA,nbIdFewNA)
            if(fast){
                resultKml <- kml3dFast(traj=traj,clusterAffectation=clustersInit)
            }else{
                if(toPlot%in%c("both","traj")){screen(listScreen[1])}else{}
                resultKml <- kml3dSlow(traj=traj,clusterAffectation=clustersInit,toPlot=toPlot,parAlgo=parAlgo)
            }

            ## A priori, une partition avec un seul cluster peut maintenant exister...
            object["add"] <- resultKml

            assign(nameObject,object,envir=parent.frame())
            if(saveCld%%parAlgo['saveFreq']==0){
                save(list=nameObject,file=paste(nameObject,".Rdata",sep=""))
                cat("S\n",saveCld," ",sep="")
            }else{
                cat("*")
            }
            if(saveCld%%100==0){cat("\n")}else{}

            if(toPlot=="both"){
                screen(listScreen[2])
                plotCriterion(as(object,"ListPartition"),nbCriterion=parAlgo['nbCriterion'])
            }else{
                if(toPlot=="criterion"){
                    plotCriterion(as(object,"ListPartition"),nbCriterion=parAlgo['nbCriterion'])
                }else{}
            }
        }
    }

    cat("\n")
    ordered(object)

    if(saveCld<Inf){
        save(list=nameObject,file=paste(nameObject,".Rdata",sep=""))
        cat("S\n")
    }else{
        cat("\n")
    }
    ## La fenetre graphique est fermée grace a 'on.exit' défini en début de fonction
    if(toPlot=="both"){
        screen(listScreen[2])
        plotCriterion(as(object,"ListPartition"),nbCriterion=parAlgo['nbCriterion'])
    }else{
        if(toPlot=="criterion"){
            plotCriterion(as(object,"ListPartition"),nbCriterion=parAlgo['nbCriterion'])
        }else{}
    }
    if(parAlgo["scale"]){restoreRealData(object)}else{}
    assign(nameObject,object,envir=parent.frame())
    return(invisible())
}

ClusterLongData3d_exportPartition3d <- function(object,nbClusters,rank,nameObject,typeGraph="bmp",parTraj=parTRAJ(),parMean=parMEAN()){
    #                           parWin=windowsCut(1)){
#    col="clusters",type="l",
#    col.mean="clusters",type.mean="b",main="",cex=1,
#    pch.mean="letters",pch.time=NA,...#,legends=TRUE,...
#){
    part <- object[paste('c',nbClusters,sep="")][[rank]]

    dataFrame <- data.frame(idAll=object["idAll"],clusters=NA)
    dataFrame$clusters[dataFrame$id%in%object['idFewNA']] <- part["clusters"]
    write.csv2(dataFrame,file=paste(nameObject,"-Clusters.csv",sep=""),row.names=FALSE)

    detail <- c(part["nbClusters"],part["percentEachCluster"],part["criterionValues"][-(length(CRITERION_NAMES)+1)],
                part["algorithm"],part["convergenceTime"])
    names(detail) <-  c("nbClusters",paste("percent",LETTERS[1:part["nbClusters"]]),CRITERION_NAMES,
                       "algorithmUsed","convergenceTime")
    write.csv2(detail,file=paste(nameObject,"-Details.csv",sep=""),row.names=TRUE)

    trajMean <- data.frame(calculTrajMean3d(object['traj'],part['clustersAsInteger']))
    write.csv2(trajMean,file=paste(nameObject,"-TrajMean.csv",sep=""),row.names=TRUE)

    eval(parse(text=paste(typeGraph,"(filename='",nameObject,"-Traj.",typeGraph,"')",sep="")))
      plotTraj(as(object,"LongData3d"),part,parTraj=parTraj,parMean=parMean)
    dev.off()
    #lty=lty,lty.mean=lty.mean,pch=pch,pch.mean=pch.mean,pch.time=pch.time,
    #xlab=xlab,ylab=ylab,ylim=ylim,cex.mean=cex.mean,legends=legends,sizeMin=sizeMin,...)
    #    savePlot(filename=paste(nameObject,"-Traj",sep=""),type=typeGraph)
    return(invisible())
}
setMethod("exportPartition",signature=c("ClusterLongData3d","numeric"),ClusterLongData3d_exportPartition3d)



#setMethod("choice",signature=c("ClusterLongData3d"),ClusterLongData_choice)



cat("\n-------------------------------------------------------------------
------------------------------ kml3d ------------------------------
------------------------------- Fin -------------------------------
-------------------------------------------------------------------\n")
