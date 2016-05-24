## clusterization est une partition associé a une longData, ou une clusterizLongData.
### cet objet ne devrait pouvoir exister que dans un cld


cat("
   ####################################################################
  ######################## Class ClustLongData3d #####################
 ############################## Creation ############################
####################################################################\n")

ClusterLongData3d_validity <- function(object){
    validObject(as(object,"LongData3d"))
    validObject(as(object,"ListPartition"))
    return(TRUE)
}

cat("### Definition ###\n")
# id       : identifiant of the individual (or lines).
# time     : real time
# varNames : nom of the variable (single now, several in the futur)
# value    : array of the trajectories. Dim 1 is individual, 2 is time, 3 is variable(s)
setClass(
    Class="ClusterLongData3d",
    contains=c("LongData3d","ListPartition"),
    validity=ClusterLongData3d_validity
)


### Data.frame ou array en 3D
clusterLongData3d <- function(traj,idAll,time,timeInData,varNames,maxNA){
    if(missing(traj)){
        return(new("ClusterLongData3d"))
    }else{}
    ## First part : set all the parameters
    if(is.data.frame(traj)){
        if(missing(idAll)){
            idAll <- traj[,1]
        }else{}
        matr <- as.matrix(traj[,sort(na.omit(unlist(timeInData)))])
        lengthTime <- length(timeInData[[1]])
        nbVar <- length(timeInData)
        traj <- array(matr[,rank(unlist(timeInData),na.last="keep")],c(nrow(traj),lengthTime,nbVar))
    }else{
        if(is.array(traj)){
            if(missing(idAll)){
                idAll <- paste("i",1:nrow(traj),sep="")
            }else{}
            if(!missing(timeInData)){
                traj <- traj[,timeInData,,drop=FALSE]
            }else{}
            lengthTime <- dim(traj)[2]
            nbVar <- dim(traj)[3]
        }else{
            stop("[ClusterLongData3d:constructor]: 'traj' should be either a data.frame or an array")
        }
    }
    if(missing(maxNA)){maxNA <- lengthTime-2}else{}
    if(length(maxNA)==1){maxNA <- rep(maxNA,nbVar)}else{}
    if(missing(varNames)){
        if(!missing(timeInData)){
            if(!is.null(names(timeInData))){
                varNames <- names(timeInData)
            }else{
                varNames <- paste("V",1:nbVar,sep="")
            }
        }else{
            varNames <- paste("V",1:nbVar,sep="")
        }
    }else{}
    if(missing(time)){time <- 1:lengthTime}else{}

    ## Second part : all the arguments are non-missing, the object can be build.

    keepId <- apply(t(apply(traj,c(1,3),function(x){sum(is.na(x))}))<=maxNA,2,all)

    traj <- traj[keepId,,,drop=FALSE]
    idFewNA <- idAll[keepId]
    dimnames(traj) <- list(idFewNA,paste("t",time,sep=""),varNames)
    reverse <- matrix(c(0,1),2,length(varNames),dimnames=list(c("mean","sd"),varNames))
    return(new("ClusterLongData3d",
        idAll=as.character(idAll),
        idFewNA=as.character(idFewNA),
        time=time,
        varNames=varNames,
        traj=traj,
        dimTraj=dim(traj),
        maxNA=maxNA,
        reverse=reverse)
    )
}

cld3d <- clusterLongData3d


ClusterLongData3d_show <- function(object){
    cat("   ~~~ Class: ClusterLongData3d ~~~")
    cat("\n      ~ Sub-Class: LongData3d ~ ")
    LongData3d_show(as(object,"LongData3d"))
    cat("\n    ~ Sub-Class: ListPartition ~ ")
    ListPartition_show(as(object,"ListPartition"))
}
setMethod("show","ClusterLongData3d",ClusterLongData3d_show)




cat("
  ###################################################################
 ############################ Get & Set ############################
###################################################################\n")

cat("### Getteur ###\n")
ClusterLongData3d_get <- function (x, i, j="missing", ..., drop = TRUE){
      if (is.numeric(i)) {
        stop("[ClusterLongData3d:getteur]: to get a clusters list, use ['ci']")
      }else{}
      if (i %in% c("criterionValues", "criterionValuesAsMatrix")){
        j <- x['criterionActif']
      }else{}
      if (i %in% c(CRITERION_NAMES, "criterionActif", CLUSTER_NAMES,
                   "criterionValues", "criterionValuesAsMatrix", "sorted",
                   "initializationMethod")) {
        x <- as(x, "ListPartition")
      }else{
        x <- as(x, "LongData3d")
      }
      return(x[i, , drop=drop])
}


setMethod(
  "[",
  signature=signature(x="ClusterLongData3d", i="ANY", j="missing",drop="ANY"),
  definition=ClusterLongData3d_get
)


cat("### Setteur ###\n")
### Héritage direct de ListPartition puisque set n'est pas défini pour LongData
### ATTENTION !!! Normalement, il faudrait vérifier que la partition est de la BONNE TAILLE !!!


ClusterLongData3d_set <- function (x, i, j="missing", ..., value){
    if (i == "add") {
      if (length(value["clusters"]) != length(x["idFewNA"])) {
        stop("[ClusterLongData3d:set] the lenght of the Partition should be the same than 'idFewNA'")
      }else{}
    }
    callNextMethod(x, i, j=j, ..., value = value)
}

setMethod(
  f="[<-",
  signature=signature(x="ClusterLongData3d", i="character", j="missing",value="missing"),
  definition=ClusterLongData3d_set
)


cat("
  #####################################################################
 ############################### plot3d ##############################
#####################################################################")

### On a un cld et un num, on plot le longData et la Partition qui va avec.
ClusterLongData3d_num_plot3d <- function(x,y,varY=1,varZ=2,parTraj=parTRAJ(),parMean=parMEAN(),...){
    if(length(y)==1){y<-c(y,1)}else{}
    yPartition <- x[paste('c',y[1],sep="")][[y[2]]]
    plotTrajMeans3d(x=as(x,"LongData3d"),y=yPartition,varY=varY,varZ=varZ,parTraj=parTraj,parMean=parMean,...)
    return(invisible())
}
setMethod("plot3d",signature=c("ClusterLongData3d","numeric"),ClusterLongData3d_num_plot3d)


### Si y est manquant :
###  - soit il est calculable et on le calcul puis on appelle plot.ClusterLongData3d
###  - soit il n'est pas calculable et on appelle plot.LongData3d.num
ClusterLongData3d_missingY_plot3d <- function(x,y,varY=1,varZ=2,parTraj=parTRAJ(),parMean=parMEAN(),...){
    plotTrajMeans3d(x=as(x,"LongData3d"),varY=varY,varZ=varZ,parTraj=parTraj,...)
    return(invisible())
}
setMethod("plot3d",signature=c("ClusterLongData3d","missing"),ClusterLongData3d_missingY_plot3d)

#setMethod("plot3d",signature=c("ClusterLongData3d","Partition"),function(x,y,...){plotTraj3d(x,y,...)})




### On a un cld et un num, on plot le longData et la Partition qui va avec.
ClusterLongData3d_num_plotTraj3d <- function(x,y,varY=1,varZ=2,parTraj=parTRAJ(col="clusters"),parMean=parMEAN(type="n"),...){
    if(length(y)==1){y<-c(y,1)}else{}
    yPartition <- x[paste('c',y[1],sep="")][[y[2]]]
    plotTrajMeans3d(x=as(x,"LongData3d"),y=yPartition,varY=varY,varZ=varZ,parTraj=parTraj,parMean=parMean,...)
    return(invisible())
}
setMethod("plotTraj3d",signature=c("ClusterLongData3d","numeric"),ClusterLongData3d_num_plotTraj3d)

ClusterLongData3d_num_plotMeans3d <- function(x,y,varY=1,varZ=2,parTraj=parTRAJ(type="n"),parMean=parMEAN(),...){
    if(length(y)==1){y<-c(y,1)}else{}
    yPartition <- x[paste('c',y[1],sep="")][[y[2]]]
    plotTrajMeans3d(x=as(x,"LongData3d"),y=yPartition,varY=varY,varZ=varZ,parTraj=parTRAJ(type="n"),parMean=parMean,...)
    return(invisible())
}
setMethod("plotMeans3d",signature=c("ClusterLongData3d","numeric"),ClusterLongData3d_num_plotMeans3d)






cat("
  #####################################################################
 ############################## plot3dPdf ############################
#####################################################################")

ClusterLongData3d_num_plot3dPdf <- function(x,y,varY=1,varZ=2){
    if(length(y)==1){y<-c(y,1)}else{}
    yPartition <- x[paste('c',y[1],sep="")][[y[2]]]
    return(plot3dPdf(x=as(x,"LongData3d"),y=yPartition,varY=varY,varZ=varZ))
}
setMethod("plot3dPdf",signature=c("ClusterLongData3d","numeric"),ClusterLongData3d_num_plot3dPdf)



cat("
  #####################################################################
 ################################ plot ###############################
#####################################################################")

### On a un cld3d et un num, on plot le longData3d et la Partition qui va avec. On ne ferme pas le screen.
ClusterLongData3d_plotTrajMeans <- function(x,y=NA,parTraj=parTRAJ(),parMean=parMEAN(),xlab="Times",ylab=x["varNames"],addLegend=TRUE,adjustLegend=-0.12,...){
    ## ############################# Preparation ############################# ##
    nbVar <- x['nbVar']
    nbTime <- x['nbTime']
    traj <- x['traj']
    parWin <- windowsCut(x['nbVar'],addLegend=addLegend)

    ## Gestion de la partition

     if(is.numeric(y)){
        if(length(y)==1){y<-c(y,1)}else{}
        y <- x[paste('c',y[1],sep="")][[y[2]]]
    }else{}

    ## Calcul du layout
    listScreen<-split.screen(parWin['screenMatrix'])

    for (i in 1:nbVar){
        screen(listScreen[i])
        par(mar=c(4,4,2,2))
        plotTrajMeans(longDataFrom3d(x,i),y,parTraj=parTraj,parMean=parMean,xlab=xlab,ylab=ylab[i],...)

    }
    if(addLegend){
        screen(listScreen[length(listScreen)],FALSE)
        part <- factor(y['clusters'],levels=LETTERS[1:y['nbClusters']])
        plotLegend(as.numeric(table(part)/length(part)*100),parMean=parMean,adjustLegend=adjustLegend)
    }

    return(listScreen)
}


##setMethod("plot",signature=c("ClusterLongData3d","missing"),.clusterLongData3d.plot)
ClusterLongData3d_plot <- function(x,y=NA,parTraj=parTRAJ(),parMean=parMEAN(),addLegend=TRUE,adjustLegend=-0.05,
    toPlot="both",nbCriterion=1000,xlab="Times",ylab=x["varNames"],closeScreenTraj=TRUE,...){
    if(any(is.na(y))){
       toPlot <- "traj"
       addLegend <- FALSE
    }else{}
    switch(EXPR=toPlot,
           "both"={
               listScreen <- split.screen(matrix(c(0,0.3,0.3,1,0,0,1,1),2))
               screen(listScreen[2])
               ClusterLongData3d_plotTrajMeans(x,y,parTraj=parTraj,parMean=parMean,addLegend=addLegend,adjustLegend=adjustLegend,xlab=xlab,ylab=ylab,...)
               screen(listScreen[1])
               plotCriterion(as(x,"ListPartition"),criterion=x["criterionActif"],nbCriterion=nbCriterion)

           },
           "traj"={
               ClusterLongData3d_plotTrajMeans(x,y=y,parTraj=parTraj,parMean=parMean,addLegend=addLegend,adjustLegend=adjustLegend,xlab=xlab,ylab=ylab,...)
           },
           "criterion"={
               plotCriterion(as(x,"ListPartition"),criterion=x['criterionActif'],nbCriterion=nbCriterion)
           }


     )
    if(closeScreenTraj){
        close.screen(all.screens = TRUE)
        return(invisible())
    }else{
        return(listScreen)
    }
}




ClusterLongData3d_missing_plot <- function(x,y,parTraj=parTRAJ(),parMean=parMEAN(),toPlot="both",nbCriterion=1000,xlab="Times",ylab=x["varNames"],...){
     ClusterLongData3d_plot(x=x,y=NA,parTraj=parTraj,parMean=parMean,toPlot=toPlot,nbCriterion=nbCriterion,xlab=xlab,ylab=ylab,closeScreenTraj=TRUE,...)
}

setMethod("plot",signature=c("ClusterLongData3d","numeric"),ClusterLongData3d_plot)
setMethod("plot",signature=c("ClusterLongData3d","missing"),ClusterLongData3d_missing_plot)


ClusterLongData3d_plotTraj <- function(x,y,xlab="Times",ylab=x["varNames"],...){
   plot(x=x,y=y,parMean=parMEAN(type="n"),toPlot="traj",parTraj=parTRAJ(col="clusters"),xlab=xlab,ylab=ylab,...)
}
setMethod("plotTraj",signature=c("ClusterLongData3d","numeric"),ClusterLongData3d_plotTraj)


ClusterLongData3d_plotMeans <- function(x,y,...){
   plot(x=x,y=y,toPlot="traj",parTraj=parTRAJ(type="n"),...)
}
setMethod("plotMeans",signature=c("ClusterLongData3d","numeric"),ClusterLongData3d_plotMeans)





gald3d <- generateArtificialLongData3d <- function(
    nbEachClusters=50,time=0:10,varNames=c("V","T"),
    meanTrajectories=list(function(t){c(0,0)},function(t){c(10,10)},function(t){c(10-t,10-t)}),
    personalVariation=function(t){c(rnorm(1,0,2),rnorm(1,0,2))},
    residualVariation=function(t){c(rnorm(1,0,2),rnorm(1,0,2))},
    decimal=2,percentOfMissing=0#,clusterLongData=TRUE
){
    nbClusters <- length(meanTrajectories)
    if(length(nbEachClusters)==1){nbEachClusters <- rep(nbEachClusters,nbClusters)}else{}
    if(is.numeric(personalVariation)){eval(parse(text=paste("personalVariation <- function(t){c(rnorm(1,0,",personalVariation,"),rnorm(1,0,",personalVariation,"))}",sep="")))}else{}
    if(length(personalVariation)==1){personalVariation <- rep(list(personalVariation),nbClusters)}else{}
    if(is.numeric(residualVariation)){eval(parse(text=paste("residualVariation <- function(t){c(rnorm(1,0,",residualVariation,"),rnorm(1,0,",residualVariation,"))}",sep="")))}else{}
    if(length(residualVariation)==1){residualVariation <- rep(list(residualVariation),nbClusters)}else{}
    if(length(percentOfMissing)==1){percentOfMissing <- rep(percentOfMissing,nbClusters)}else{}
    nbTime <- length(time)
    nbVar <- length(varNames)
    idAll <- paste("i",1:(sum(nbEachClusters)),sep="")
    indivInCluster <- rep(1:nbClusters,times=nbEachClusters)

    traj <- array(NA,dim=c(sum(nbEachClusters),nbTime,nbVar),dimnames=list(idAll,paste("t",time,sep=""),varNames))
    for (iIndiv in 1:nrow(traj)){
        traj[iIndiv,,] <- t(sapply(time,meanTrajectories[[indivInCluster[iIndiv]]])+personalVariation[[indivInCluster[iIndiv]]](0)+sapply(time,residualVariation[[indivInCluster[iIndiv]]]))
    }
    traj <- round(traj,digits=decimal)

    for (iCluster in 1:nbClusters){
        nbVal <- nbTime*nbEachClusters[iCluster]
        while(sum(is.na(traj[indivInCluster==iCluster,,]))/nbVal < percentOfMissing[iCluster]){
            randL <- floor(runif(1,cumsum(c(0,nbEachClusters))[iCluster]+1,cumsum(nbEachClusters)[iCluster]+1))
            randC <- floor(runif(1,1,nbTime+1))
            randV <- floor(runif(1,1,nbVar+1))
            if(sum(!is.na(traj[randL,,randV]))>1){traj[randL,randC,randV]<-NA}else{}
        }
    }
#    if(clusterLongData){return(as.clusterLongData(traj,idAll=id,time=time,varNames=varNames))}else{
    return(clusterLongData3d(traj,idAll=idAll,time=time,varNames=varNames))
}




cat("\n--------------------------------------------------------------------
------------------------ Class ClustLongData3d ---------------------
--------------------------------- Fin ------------------------------
--------------------------------------------------------------------\n")
