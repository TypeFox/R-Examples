cat("\n####################################################################
################################ kml ###############################
############################# Creation #############################
####################################################################\n")



### On suppose que si un centre est NA, il est en dernière ligne de clustersCenter
calculTrajFuzzyMean <- function(traj,fuzzyClust){
    nbTime <- ncol(traj)
    nbClust <- ncol(fuzzyClust)
    trajMean <- matrix(NA,nbClust,nbTime)
    for(k in 1:nbClust){
        trajMean[k,] <- apply(traj,2,function(x){weighted.mean(x,fuzzyClust[,k])})
    }
    return(trajMean)
}



affectFuzzyIndiv <- function(traj,clustersCenter,fuzzyfier=1.25){
    nbIndiv <- nrow(traj)
    nbClust <- nrow(clustersCenter)
    clusterAffectation <- matrix(0,nbIndiv,nbClust)

    for(i in 1:nbIndiv){
        distCenter <- numeric(nbClust)
        distCenter <- apply(clustersCenter,1,function(x)dist(rbind(traj[i,],x))^(2/(fuzzyfier-1)))
        if(min(distCenter)==0){
            clusterAffectation[i,which.min(distCenter)] <- 1
            ## Les autres sont a zéro par défaut
        }else{
            clusterAffectation[i,] <- 1/(distCenter*sum(1/distCenter))
        }
    }
    return(clusterAffectation)
}


fuzzyKmlSlow <- function(traj,clusterAffectation,toPlot="traj",fuzzyfier=1.25,parAlgo=parALGO()){
    nbId <- nrow(traj)
    nbClusters <- max(clusterAffectation,na.rm=TRUE)

    clusterAppartenance <- matrix(0,nbId,max(clusterAffectation,na.rm=TRUE))
    clusterAppartenance[na.omit(matrix(c(1:nbId,clusterAffectation),nbId))] <- 1
    exClusterAffectation <- exExClusterAffectation <- exExExClusterAffectation <- NA
    if(toPlot%in%c("traj","both")){
        plotTrajMeans(longData(traj),partition(clusterAffectation))
    }else{}
    for(iterations in 1:parAlgo['maxIt']){
        clustersCenter <- calculTrajFuzzyMean(traj=traj,fuzzyClust=clusterAppartenance)
        clusterAppartenance <- affectFuzzyIndiv(traj=traj,clustersCenter=clustersCenter,fuzzyfier=fuzzyfier)
        clusterAffectation <- apply(clusterAppartenance,1,which.max)
        if(toPlot%in%c("traj","both")){
#            plot(longData(traj),partition(clusterAffectation))
            matplot(t(traj),type="l",col=1)
            matlines(t(clustersCenter),col=1:nbClusters+1,type="l",lwd=7,lty=1)
        }else{}
        if(identical(clusterAffectation,exExExClusterAffectation)){
            return(clusterAppartenance)
        }else{
            exExExClusterAffectation <- exExClusterAffectation
            exExClusterAffectation <- exClusterAffectation
            exClusterAffectation <- clusterAffectation
        }
    }
    return(clusterAppartenance)
}


### ATTENTION : calculCenterGeneralized travaillait avec traj + Partition ;
### Maintenant, il travaille avec traj + part.
#calculCenterGeneralized <- calculTrajMean <- function(traj,part,centerMethod=meanNA){
#    trajMean <- apply(traj, c(2,3), tapply, part, centerMethod)
#    return(trajMean)
#}

### ATTENTION : Ne fonctionne pas avec les partitions à un seul clusters
calculTrajMean <- function(traj,clust,centerMethod=function(x){mean(x,na.rm=TRUE)}){
    trajMean <- apply(traj, 2, tapply, clust, centerMethod)
    return(trajMean)
}


### On suppose que si un centre est NA, il est en dernière ligne de clustersCenter
affectIndiv <- function(traj,clustersCenter,distance=function(x,y){dist(rbind(x,y))}){
#    if (distance %in% METHODS){distanceFun <- ,method=distance))}}else{distanceFun <- distance}
    nbId <- nrow(traj)
    clusterAffectation <- rep(1,nbId)
    distActuel <- apply(traj,1,function(x){distance(x,clustersCenter[1,])})
    ##   print(distActuel)
    for(iNbClusters in 2:nrow(clustersCenter)){
        distToMean <- apply(traj,1,function(x){distance(x,clustersCenter[iNbClusters,])})
 #       print(distToMean)
        cond <- distToMean<distActuel
        cond[is.na(cond)] <- FALSE # Car si cond==NA, c'est que distToMean==NA et donc on ne change pas l'affectation.
        clusterAffectation <- ifelse(cond,rep(iNbClusters,nbId),clusterAffectation)
        distActuel <- ifelse(distToMean<distActuel,distToMean,distActuel)
    }
    return(clusterAffectation)
}


calculTrajMeanC  <- function(traj,clust){
    # print(as.integer(xPart["clusters"]))
    # L'initialisation de trajMean sous R a 0 est indispensable.
    nbClusters <- max(clust,na.rm=TRUE)
    trajMean <- matrix(0,nbClusters,ncol(traj))
    result <- .C("calculMean",traj=as.double(t(traj)),nbInd=as.integer(nrow(traj)),nbTime=as.integer(ncol(traj)),
        clusterAffectation=as.integer(clust),nbCluster=as.integer(nbClusters),trajMean=as.numeric(t(trajMean)),NAOK=TRUE,PACKAGE="kml")$trajMean
    return(matrix(result,nbClusters,byrow=TRUE))
}


affectIndivC <- function(traj,clustersCenter){
    part <- rep(0,nrow(traj))
    result <- .C("affecteIndiv",traj=as.double(t(traj)),nbInd=as.integer(nrow(traj)),nbTime=as.integer(ncol(traj)),
        trajMean=as.numeric(t(clustersCenter)),nbCluster=as.integer(nrow(clustersCenter)),clusterAffectation=as.integer(part),NAOK=TRUE,PACKAGE="kml"
    )
#    print(result)
    return(result$clusterAffectation)
}


kmlSlow <- function(traj,clusterAffectation,toPlot="traj",parAlgo=parALGO()){
#    if (distance %in% METHODS){distanceFun <- function(x,y){return(dist(t(cbind(x,y)),method=distance))}}else{distanceFun <- distance}
 #   print(distanceFun)
    longDatTraj <- longData(traj,maxNA=ncol(traj)-1)
    kmlCenterMethod <- parAlgo['centerMethod']
    kmlDistance <- parAlgo['distance']

    exClusterAffectation <- NA

    if(toPlot==c("traj")){
        ClusterLongData_plotTrajMeans(longDatTraj,partition(clusterAffectation),addLegend=TRUE,xlab="Times",ylab="V")
    }else{}
    if(toPlot==c("both")){
        screen(2)
        ClusterLongData_plotTrajMeans(longDatTraj,partition(clusterAffectation),addLegend=TRUE,xlab="Times",ylab="V")
    }else{}

    for(iterations in 1:parAlgo['maxIt']){
        clustersCenter <- calculTrajMean(traj=traj,clust=clusterAffectation,centerMethod=kmlCenterMethod)
        clusterAffectation <- affectIndiv(traj=traj,clustersCenter=clustersCenter,distance=kmlDistance)

        if(toPlot==c("traj")){
            ClusterLongData_plotTrajMeans(longDatTraj,partition(clusterAffectation),addLegend=TRUE,xlab="Times",ylab="V")
        }else{}
        if(toPlot==c("both")){
            screen(2)
            ClusterLongData_plotTrajMeans(longDatTraj,partition(clusterAffectation),addLegend=TRUE,xlab="Times",ylab="V")
        }else{}

        if(identical(clusterAffectation,exClusterAffectation)){
            clusterAffectation <- partition(clusterAffectation,longDatTraj,
               details=c(convergenceTime=as.character(iterations),algorithm="kmeans, slow (R)",multiplicity="1"))
            return(clusterAffectation)
        }else{
            exClusterAffectation <- clusterAffectation
        }
    }
    return(partition(clusterAffectation,longDatTraj,details=c(convergenceTime=as.character(Inf),algorithm="kmeans, slow (R)",multiplicity="1")))
}


#traj <- ld3['traj']
#clusterAffectation <- partitionInitialise(3,8,method="randomAll")

kmlFast <- function(traj,clusterAffectation){
    resultKml <- .C("kml1",as.double(t(traj)),nbInd=as.integer(nrow(traj)),nbTime=as.integer(ncol(traj)),
                    nbCluster=as.integer(max(clusterAffectation,na.rm=TRUE)),maxIt=as.integer(200),
                    clusterAffectation1=as.integer(clusterAffectation),
                    convergenceTime=as.integer(1),
                    NAOK=TRUE,PACKAGE="kml")[c(6,7)]
    return(partition(resultKml[[1]],longData(traj,maxNA=ncol(traj)-1),details=c(convergenceTime=as.character(resultKml[[2]]),algorithm="kmeans, fast (C)",multiplicity="1")))
}


expandStartingCond <- function(startingCond,nbRedrawing,methodUsed){
    startingSeq <- character()
    if(length(startingCond)==1){
        if(startingCond%in%c("all","nearlyAll")){
            if(startingCond=="all"&!"maxDist"%in%methodUsed){startingSeq <- c(startingSeq,"maxDist")}else{}
            if(startingCond%in%c("all","nearlyAll")&!"kmeans-"%in%methodUsed){startingSeq <- c(startingSeq,"kmeans-")}else{}
            startingSeq <- c(startingSeq,rep(c("kmeans--","randomK"),nbRedrawing))
        }else{
            startingSeq <- rep(startingCond,nbRedrawing)
        }
    }else{
        startingSeq <- rep(startingCond,nbRedrawing)
    }
    return(startingSeq[1:nbRedrawing])
}

### Si on ne doit faire qu'un seul plot, on coute tout de même en deux
###   (deux régions qui ont la même taille, tout le dessin)
### Ca permet de definir screen(2) pour les criterions dans tous les cas
###   et screen(1) pour les traj dans tous les cas.

cutScreen <- function(toPlot){
    return(switch(EXPR=toPlot,
          "both"={split.screen(matrix(c(0.3,0,1,0.3,0,0,1,1),2))},
          "traj"={split.screen(matrix(c(0,0,1,1,0,0,1,1),2))},
          "criterion"={split.screen(matrix(c(0,0,1,1,0,0,1,1),2))},
          "none"=NULL)
    )
}



#exCutScreen <- function(toPlot){
#        if(toPlot=="both"){
#        listScreen <- split.screen(matrix(c(0.3,0,1,0.3,0,0,1,1),2))
#    }else{
#        if(
#        listScreen <- split.screen(c(1,1))
#    }
#    return(listScreen)
#}

fastOrSlow <- function(toPlot,parAlgo){
    if((toPlot%in%c("criterion","none")) & (parAlgo['distanceName']=="euclidean") & identical(parAlgo["centerMethod"],parALGO()["centerMethod"])){
        cat(" ~ Fast KmL ~\n")
        fast <- TRUE
    }else{
        cat(" ~ Slow KmL ~\n")
        fast <- FALSE
    }
    return(fast)
}




kml <- function(object,nbClusters=2:6,nbRedrawing=20,toPlot="none",parAlgo=parALGO()){
    if(class(object)=="ClusterLongData3d"){
        stop("[kml]: kml is for longitudinal data (object 'ClusterLongData').
For joint longitudinal data (object of class 'ClusterLongData3d'), use kml3d")
    }else{}

    nameObject<-deparse(substitute(object))
    on.exit(if(toPlot!="none"){close.screen(all.screens=TRUE)}else{})

    nbIdFewNA <- object["nbIdFewNA"]
    convergenceTime <- 0
    traj <- object["traj"]
    nbTime <- length(object["time"])
    saveCld <-0

    ################
    ## listScreen[1] (à droite) est pour les traj.
    if(toPlot!="none"){
        plot(object,closeScreenTraj=(toPlot!="both"))
    }else{}

    ################
    ## Starting conditions
    startingCond <- expandStartingCond(parAlgo['startingCond'],nbRedrawing,object["initializationMethod"])
    object["initializationMethod"] <- unique(c(object["initializationMethod"],startingCond))

    ################
    ## Fast or Slow, according to parAlgo to toPlot
    fast <- fastOrSlow(toPlot,parAlgo)

    for(iRedraw in 1:nbRedrawing){
        for(iNbClusters in nbClusters){
            saveCld <- saveCld+1
            clustersInit <- initializePartition(nbClusters=iNbClusters,lengthPart=nbIdFewNA,method=startingCond[iRedraw],data=traj)
            clust <- rep(NA,nbIdFewNA)
            if(fast){
                resultKml <- kmlFast(traj=traj,clusterAffectation=clustersInit)
            }else{
                resultKml <- kmlSlow(traj=traj,clusterAffectation=clustersInit,toPlot=toPlot,parAlgo=parAlgo)
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
            if(toPlot=="both"){
                screen(1)
                plotCriterion(as(object,"ListPartition"),nbCriterion=parAlgo['nbCriterion'])
            }else{
                if(toPlot=="criterion"){
                    plotCriterion(as(object,"ListPartition"),nbCriterion=parAlgo['nbCriterion'])
                }else{}
            }
        }
    }
    ## La fenetre graphique est fermée grace a 'on.exit' défini en début de fonction
    ordered(object)
    if(parAlgo["saveFreq"]<Inf){
        save(list=nameObject,file=paste(nameObject,".Rdata",sep=""))
        cat("S\n")
    }else{
        cat("\n")
    }

    if(toPlot=="both"){
        screen(1)
        plotCriterion(as(object,"ListPartition"),nbCriterion=parAlgo['nbCriterion'])
    }else{
        if(toPlot=="criterion"){
            plotCriterion(as(object,"ListPartition"),nbCriterion=parAlgo['nbCriterion'])
        }else{}
    }
    assign(nameObject,object,envir=parent.frame())
    return(invisible())
}




exportPartition <- function(object,nbClusters,rank,nameObject,typeGraph="bmp",parTraj=parTRAJ(),parMean=parMEAN()){
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

    trajMean <- data.frame(calculTrajMeanC(object['traj'],part['clustersAsInteger']))
    write.csv2(trajMean,file=paste(nameObject,"-TrajMean.csv",sep=""),row.names=TRUE)

    eval(parse(text=paste(typeGraph,"(filename='",nameObject,"-Traj.",typeGraph,"')",sep="")))
    plotTrajMeans(as(object,"LongData"),part,parTraj=parTraj,parMean=parMean)
    dev.off()
        #lty=lty,lty.mean=lty.mean,pch=pch,pch.mean=pch.mean,pch.time=pch.time,
        #xlab=xlab,ylab=ylab,ylim=ylim,cex.mean=cex.mean,legends=legends,sizeMin=sizeMin,...)
#    savePlot(filename=paste(nameObject,"-Traj",sep=""),type=typeGraph)
    return(invisible())
}

#setMethod("exportPartition",signature=c("ClusterLongData","numeric"),.exportPartition)

#exportPartition(cld4,3,1,"testPart")




choiceChangeParam <- function(paramChoice){
    xy <- paramChoice['xy']

    texte <- paste("     ~ Choice : menu ~
 - 'Arrow' : change partition
 - 'Space' : select/unselect a partition
 -    e    : change the display (",paramChoice['toPlot'],")
 -   d/s   : change actif criterion (",CRITERION_NAMES[paramChoice['critRank']],")
 -    c    : sort according to the actif criterion (",paramChoice['critSorted'][paramChoice['critRank']],"),
 -    r    : change the trajectories's style (type=",CHOICE_STYLE[['typeTraj']][paramChoice['styleTrajRank']],
                                            "; col=",CHOICE_STYLE[['colTraj']][paramChoice['styleTrajRank']],")
 -    f    : change the means trajectories' style (type=",CHOICE_STYLE[['typeMean']][paramChoice['styleMeanRank']],
                                    "; col=",CHOICE_STYLE[['colMean']][paramChoice['styleMeanRank']],
                                    "; pch=",CHOICE_STYLE[['pchMean']][paramChoice['styleMeanRank']],")
 -   g/t   : change the symbol size (",paramChoice['cex'],")
 -   y/h   : change the number of symbols (freq=1/",1+paramChoice['pchPeriod'],")
 -   j/u   : legend down/up (",paramChoice['yLegend'],")

     ~ 'Return' or 'm' when its done ~\n",sep="")

    choix <- getGraphicsEvent(texte,onKeybd=function(key){return(key)})
    cat("Choix : ",choix," class :",class(choix)," length :",length(choix),"\n")
    switch(EXPR=choix,
           "Up"    = {
               if(xy[1]>1){
                   paramChoice['toDo'] <- "xy"
                   xy[2]<-1
                   xy[1]<-xy[1]-1
                   paramChoice['xy']<-xy
#                       if(is.na(critMatrix[y[1],1])){
 #                          y[1] <- y[1]+1-which.min(is.na(critMatrix[,1][y[1]:1]))
  #                         if(is.na(critMatrix[y[1],1])){
   #                            y[1] <- y[1]-1+which.min(is.na(critMatrix[,1][y[1]:52]))
    #                       }else{}
     #                  }else{}
               }else{paramChoice['toDo'] <- ""}
           },
           "Down"  = {
               if(xy[1]<nrow(paramChoice['critMatrix'])){
                   paramChoice['toDo'] <- "xy"
                   xy[2]<-1
                   xy[1]<-xy[1]+1
                   paramChoice['xy']<-xy
#                       if(is.na(critMatrix[y[1],1])){
 #                          y[1] <- y[1]-1+which.min(is.na(critMatrix[,1][y[1]:52]))
  #                         if(is.na(critMatrix[y[1],1])){
   #                            y[1] <- y[1]+1-which.min(is.na(critMatrix[,1][y[1]:1]))
    #                       }else{}
     #                  }else{}
               }else{paramChoice['toDo'] <- ""}
           },
           "Right" = {
               paramChoice['toDo'] <- "xy"
               if(xy[2]<ncol(paramChoice['critMatrix']) && !is.na(paramChoice['critMatrix'][xy[1],xy[2]+1])){
                   paramChoice['toDo'] <- "xy"
                   xy[2]<-xy[2]+1
                   paramChoice['xy']<-xy
               }else{paramChoice['toDo'] <- ""}
           },
           "Left"  = {
               paramChoice['toDo'] <- "xy"
               if(xy[2]>1){
                   paramChoice['toDo'] <- "xy"
                   xy[2]<-xy[2]-1
                   paramChoice['xy']<-xy
               }else{paramChoice['toDo'] <- ""}

           },

           "ctrl-J" = {
               paramChoice['toDo'] <- "EXIT"
           },
           "m" = {
               paramChoice['toDo'] <- "EXIT"
           },
           " "      = {
               paramChoice['toDo'] <- ""
               if(list(xy) %in% paramChoice['selectedPart']){
                   paramChoice['selectedPart'] <- paramChoice['selectedPart'][!(paramChoice['selectedPart'] %in% list(xy))]
               }else{
                   paramChoice['selectedPart'] <- c(paramChoice['selectedPart'],list(xy))
               }
           },
           "e" = {
               paramChoice['toDo'] <- ""
               paramChoice['toPlot'] <- ifelse(paramChoice['toPlot']=="both","traj",
                                               ifelse(paramChoice['toPlot']=="traj","criterion","both")
                                               )
           },
           "r" = {
               paramChoice['toDo'] <- "parTraj"
               paramChoice['styleTrajRank'] <- paramChoice['styleTrajRank']%%3+1
           },
           "f" = {
               paramChoice['toDo'] <- "parMean"
               paramChoice['styleMeanRank'] <- paramChoice['styleMeanRank']%%7+1
           },
           "t" = {
               paramChoice['toDo'] <- "parMean"
               paramChoice['cex'] <- paramChoice['cex']+0.1
           },
           "g" = {
               paramChoice['toDo'] <- "parMean"
               paramChoice['cex'] <- paramChoice['cex']-0.1
           },
           "h" = {
               paramChoice['toDo'] <- "parMean"
               paramChoice['pchPeriod'] <- ceiling(paramChoice['pchPeriod']*1.05+0.05)
#               if(paramChoice['pchPeriod']>paramChoice['nbTime']){
 #                  paramChoice['pchPeriod'] <- paramChoice['nbTime']}else{}
           },
           "y" = {
               paramChoice['toDo'] <- "parMean"
               paramChoice['pchPeriod'] <- floor(paramChoice['pchPeriod']/1.05)
               if(paramChoice['pchPeriod']<0){paramChoice['pchPeriod'] <- 0}else{}
           },
           "d" = {
               paramChoice['toDo'] <- "changeCriterion"
               paramChoice['critRank'] <- (paramChoice['critRank']%%length(CRITERION_NAMES)) + 1
           },
           "s" = {
               paramChoice['toDo'] <- "changeCriterion"
               paramChoice['critRank'] <- ((paramChoice['critRank']-2)%%length(CRITERION_NAMES)) + 1
#               object['criterionActif'] <- paramChoice['critPossible'][paramChoice['critRank']]
 #              ordered(object)
#               critMatrix <- object["criterionValues",paramChoice['critPossible'][paramChoice['critRank']]]
 #              lengthList <- max(sapply(critMatrix , length))
  #             critMatrix <- t(sapply(critMatrix , function(x) c(x,rep(NA,lengthList-length(x)))))
   #            paramChoice['critMatrix'] <- critMatrix
           },
           "c" = {
               paramChoice['toDo'] <- "order"
           },
           "j" = {
               paramChoice['toDo'] <- ""
               paramChoice["yLegend"] <- paramChoice["yLegend"]+0.01
           },
           "u" = {
               paramChoice['toDo'] <- ""
               paramChoice["yLegend"] <- paramChoice["yLegend"]-0.01
           },
           default={}

           )
    return(paramChoice)
}
#paramChoice <- parChoice()
#cleanProg(choiceChangeParam,,,2) # CHOICE_STYLE length



partPermut <- function(selectedPart,matPermut){
    onePermut <- function(xy){
        xy[2] <- which(matPermut[xy[1],]%in%xy[2])
        return(xy)
    }
    return(lapply(selectedPart,onePermut))
}

cat("### Method: 'choice' pour clusterizLongData ###\n")
#ClusterLongData_
choice <- function(object,typeGraph="bmp"){
    nameObject <- deparse(substitute(object))

   # pchTime <- object["time"]
#    pchFreq <- length(pch.time)

  #  size <- 1
 #   nbTime <- object["nbTime"]

    ## Fonction qui 'cercle' les selected
    pointCal <- function(z){points(z[2],critMatrix[z[1],z[2]],lwd=3,cex=3)}

    nbVar <- object['nbVar']
    critMatrix <- object["criterionValuesAsMatrix"]
    if(is.null(rownames(critMatrix))){critMatrix <- t(critMatrix)}else{}
    critMatrixRowName <- match(rownames(critMatrix),CLUSTER_NAMES)

    y <- as.numeric(c(which.max(critMatrix[,1]),1))
    paramTraj <- parTRAJ()
    paramMean <- parMEAN()
#    calSelected <- list()
    paramChoice <- parChoice(xy=y,nbTime=object['nbTime'],critMatrix=critMatrix,selectedPart=list(),
                             critSorted=((CRITERION_NAMES%in%object['criterionActif']) & object['sorted']))

    listScreen <- plot(object,c(critMatrixRowName[y[1]],y[2]),toPlot="both",closeScreenTraj=FALSE)
    if(paramChoice['toPlot']=="both"){points(y[2],critMatrix[y[1],y[2]],pch=19,lwd=5)}else{}
    close.screen(,TRUE)

    while(TRUE){
    #   print(paramChoice['selectedPart'])
        paramChoice <- choiceChangeParam(paramChoice)
        switch(EXPR=paramChoice['toDo'],
               "xy"={y <- paramChoice['xy']},
               "parTraj"={paramTraj <- parTRAJ(type=CHOICE_STYLE[['typeTraj']][paramChoice['styleTrajRank']],
                                               col=CHOICE_STYLE[['colTraj']][paramChoice['styleTrajRank']])
               },
               "parMean"={paramMean <- parMEAN(type=CHOICE_STYLE[['typeMean']][paramChoice['styleMeanRank']],
                                               col=CHOICE_STYLE[['colMean']][paramChoice['styleMeanRank']],
                                               pch=CHOICE_STYLE[['pchMean']][paramChoice['styleMeanRank']],
                                               pchPeriod=paramChoice['pchPeriod'],
                                               cex=paramChoice['cex'])
               },
               "changeCriterion"={
                   if(paramChoice['critRank']!=0){object['criterionActif'] <- CRITERION_NAMES[paramChoice['critRank']]}else{}
                   paramChoice['critMatrix'] <- critMatrix <- object["criterionValuesAsMatrix"]
               },
               "order"={
                   matPermut <- ordered(object)
                   paramChoice['critMatrix'] <- critMatrix <- object["criterionValuesAsMatrix"]
                   paramChoice['selectedPart'] <- partPermut(paramChoice['selectedPart'],matPermut)
                   y <- paramChoice['xy'] <- unlist(partPermut(list(paramChoice['xy']),matPermut))
               },
               "EXIT"={break;}
        )

        listScreen <- plot(object,c(critMatrixRowName[y[1]],y[2]),toPlot=paramChoice['toPlot'],parTraj=paramTraj,parMean=paramMean,closeScreenTraj=FALSE,adjustLegend=paramChoice["yLegend"])
        if(paramChoice['toPlot']=="both"){
            points(y[2],critMatrix[y[1],y[2]],pch=19,lwd=5)
            lapply(paramChoice['selectedPart'],pointCal)
        }else{}
        close.screen(,TRUE)
    }

#    nameObject <- paste(nameObject,"-C",y[1],"-",y[2],sep="")
    if(length(paramChoice['selectedPart'])!=0){
        for(iY in paramChoice['selectedPart']){
            exportPartition(object=object,nbClusters=critMatrixRowName[iY[1]],rank=iY[2],
                             nameObject=paste(nameObject,"-C",critMatrixRowName[iY[1]],"-",iY[2],sep=""),
                             typeGraph=typeGraph,parTraj=paramTraj,parMean=paramMean)
        }
        eval(parse(text=paste(typeGraph,"(filename='",nameObject,"-criterionActif.",typeGraph,"')",sep="")))
            plot(object,toPlot="criterion")
        dev.off()
        eval(parse(text=paste(typeGraph,"(filename='",nameObject,"-criterionAll.",typeGraph,"')",sep="")))
            plotAllCriterion(as(object,"ListPartition"),criterion=CRITERION_NAMES,standardized=TRUE)
        dev.off()
    }
    assign(nameObject,object,envir=parent.frame())
    return(invisible())

#    lapply(paramChoice['selectedPart'],exportSelected)
#                               col=colTrajPossible[styleTraj],type=typeTrajPossible[styleTraj],
 #                              col.mean=colMeanPossible[styleMeanTraj],type.mean=typeMeanPossible[styleMeanTraj],main="",cex=size,
  #                             pch.mean=pchMeanPossible[styleMeanTraj],pchTime=pchTime,
   #                            col.sub=colTrajPossible[styleTraj],type.sub=typeTrajPossible[styleTraj],
    #                           col.mean.sub=colMeanPossible[styleMeanTraj],type.mean.sub=typeMeanPossible[styleMeanTraj],main.sub="")
}

#setMethod("choice",signature=c("ClusterLongData"),ClusterLongData_choice)



cat("\n-------------------------------------------------------------------
------------------------------- kml -------------------------------
------------------------------- Fin -------------------------------
-------------------------------------------------------------------\n")
