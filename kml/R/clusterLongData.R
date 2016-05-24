### clusterization est une partition associé a une longData, ou une clusterizLongData.
### cet objet ne devrait pouvoir exister que dans un cld


cat("
   ####################################################################
  ######################### Class ClustLongData ######################
 ############################## Creation ############################
####################################################################\n")

ClusterLongData_validity <- function(object){
    validObject(as(object,"LongData"))
    validObject(as(object,"ListPartition"))
    return(TRUE)
}
cat("### Definition ###\n")
# id       : identifiant of the individual (or lines).
# time     : real time
# varNames : nom of the variable (single now, several in the futur)
# value    : array of the trajectories. Dim 1 is individual, 2 is time, 3 is variable(s)
setClass(
    Class="ClusterLongData",
    contains=c("LongData","ListPartition"),
    validity=ClusterLongData_validity
)


##############################
### Code copier intégralement depuis "LongData.r"

### Data.frame ou array en 2D
cld <- clusterLongData <- function(traj,idAll,time,timeInData,varNames,maxNA){

    if(missing(traj)){
        return(new("ClusterLongData"))
    }else{}

    ## First part : set all the parameters
    if(is.data.frame(traj)){
        if(missing(idAll)){
            idAll <- traj[,1]
            if(missing(timeInData)){
                timeInData <- 2:ncol(traj)
            }else{}
        }else{
            if(missing(timeInData)){
                timeInData <- 1:ncol(traj)
            }else{}
        }
        traj <- as.matrix(traj[,timeInData])
    }else{
        if(is.array(traj)){
            if(missing(idAll)){
                idAll <- paste("i",1:nrow(traj),sep="")
            }else{}
            if(missing(timeInData)){
                timeInData <- 1:ncol(traj)
            }else{}
            traj <- traj[,timeInData,drop=FALSE]
        }else{
            stop("[ClusterLongData:constructor]: 'traj' should be either a data.frame, a matrix or an array")
        }
    }
    if(missing(varNames)){varNames <- "V"}else{}
    if(missing(time)){time <- 1:ncol(traj)}else{}
    if(missing(maxNA)){maxNA <- ncol(traj)-2}else{}

    ## Second part : all the arguments are non-missing, the object can be build.

    ## X1 <- apply(traj,c(1,3),function(x){sum(is.na(x))}) compte le nombre de NA par indiv et par variable
    ## X2 <- t(X1)<=maxNA pour chaque ligne (ie chaque variable), indique TRUE si le nombre de NA est plus petit que le maxNA correspondant
    ## apply(X2,2,all) vérifie que la condition est bonne pour toutes les variables.

    keepId <- apply(t(apply(traj,1,function(x){sum(is.na(x))}))<=maxNA,2,all)

    ## Si on permet l'excusion globale, la formule est :
    ## keepId <- apply(traj,1,function(x)(sum(is.na(x))<=maxNA))
    traj <- traj[keepId,,drop=FALSE]
    idFewNA <- idAll[keepId]
    dimnames(traj) <- list(idFewNA,paste("t",time,sep=""))
    reverse <- matrix(c(0,1),2,1,dimnames=list(c("mean","sd"),varNames))
    return(new("ClusterLongData",
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


ClusterLongData_show <- function(object){
    cat("   ~~~ Class: ClusterLongData ~~~")
    cat("\n      ~ Sub-Class: LongData ~ ")
    LongData_show(as(object,"LongData"))
    cat("\n    ~ Sub-Class: ListPartition ~ ")
    ListPartition_show(as(object,"ListPartition"))
}
setMethod("show","ClusterLongData",ClusterLongData_show)





cat("### Getteur ###\n")

ClusterLongData_get <- function (x, i, j, ..., drop = TRUE){
#    .local <- function (x, i, j, drop){
    if (is.numeric(i)) {
        stop("[ClusterLongData:getteur]: to get a clusters list, use ['ci']")
    }else{}
    if (i %in% c("criterionValues", "criterionValuesAsMatrix")){
        if(missing(j)){j <- x['criterionActif']}else{}
        result <- as(x, "ListPartition")[i,j,drop=drop]
    }else{
        if (i %in% c(CRITERION_NAMES, "criterionActif", CLUSTER_NAMES,
                     "sorted","initializationMethod")) {
            result <- as(x, "ListPartition")[i,,drop=drop]
        }else{
            result <- as(x, "LongData")[i,,drop=drop]
        }
    }
    return(result)

                                        #    }.local(x, i, j, ..., drop)
}


############### Pourquoi séparer "character" de "numeric" et ne pas mettre "any" ?
#setMethod(
#  "[",
#  signature=signature(x="ClusterLongData", i="character", j="ANY",drop="ANY"),
#  definition=ClusterLongData_get
#)

#setMethod(
#  "[",
#  signature=signature(x="ClusterLongData", i="numeric", j="ANY",drop="ANY"),
#  definition=ClusterLongData_get
#)

setMethod(
  "[",
  signature=signature(x="ClusterLongData", i="ANY", j="ANY" ,drop="ANY"),
  definition=ClusterLongData_get
)


getClusters <- function(xCld,nbCluster,clusterRank=1,asInteger=FALSE){
    cluster <-  xCld["idAll"] %in% xCld["idFewNA"]
    cluster[cluster] <- xCld[paste("c",nbCluster,sep="")][[clusterRank]]["clustersAsInteger"]
    cluster[!cluster] <- NA
    if(!asInteger){cluster <- factor(LETTERS[cluster])}else{}
    return(cluster)
}


getBestPostProba <- function(xCld,nbCluster,clusterRank=1){
    bestPP <-  xCld["idAll"] %in% xCld["idFewNA"]
    bestPP[!bestPP] <- NA
    bestPP[!is.na(bestPP)] <- apply(xCld[paste("c",nbCluster,sep="")][[clusterRank]]["postProba"],1,max,na.rm=TRUE)
    return(bestPP)
}



cat("### Setteur ###\n")
### Héritage direct de ListPartition puisque set n'est pas défini pour LongData
### ATTENTION !!! Normalement, il faudrait vérifier que la partition est de la BONNE TAILLE !!!

setMethod(
  f="[<-",
  signature=signature(x="ClusterLongData", i="character", j="missing",value="missing"),
  definition=function (x, i, j="missing", ..., value){
    if (i == "add") {
      if (length(value["clusters"]) != length(x["idFewNA"])) {
        stop("[ClusterLongData:set] the lenght of the Partition should be the same than 'idFewNA'")
      }else{}
    }
    callNextMethod(x, i, j=j, ..., value = value)
  }
)

setMethod("is.na", "ClusterLongData", function(x) FALSE)



cat("\n
   ####################################################################
  ######################### Class ClustLongData ######################
 ############################### Autres #############################
####################################################################\n")



legendCol <- function(nbVar){
    if(nbVar<6){return(nbVar)}else{
        if(nbVar %in% c(6)){return(3)}else{
            if(nbVar %in% c(7,8,11,12)){return(4)}else{
                if(nbVar %in% c(9,10,13:15)){return(5)}else{
                    if(nbVar %in% c(16:18,21:24)){return(6)}else{
                        return(7)}}}}}
}


plotLegend <- function(percent,adjustLegend=-0.12,parMean=parMEAN()){
    nbClusters <- length(percent)
    parMean <- expandParLongData(parMean,nbClusters)
    percent <- paste(": ",formatC(percent,digits=3),"%",sep="")
#    par(mar=c(0,0,0,0))
 #   plot(1,axes=FALSE,type="n",xlab="",ylab="",xlim=c(-1,1),ylim=c(-1,1))
    legend("top",legend=percent,lty=1,col=parMean['col'],pch=parMean['pch'],
               ncol=legendCol(nbClusters),xpd=NA,inset=adjustLegend)
}



ClusterLongData_plotTrajMeans <- function(x,y=NA,parTraj=parTRAJ(),parMean=parMEAN(),addLegend=TRUE,adjustLegend=-0.12,...){
    ## Comme il n'y a pas de partition, 'clusters' devient 'black'
    ## Calcul du layout
## cat("\nA",screen())
    if(is.numeric(y)){
        if(length(y)==1){y<-c(y,1)}else{}
        y <- x[paste('c',y[1],sep="")][[y[2]]]
    }else{}
## cat("\nB",screen())
    plotTrajMeans(x=x,y=y,parTraj=parTraj,parMean=parMean,...)
## cat("\nC",screen())
    if(addLegend&!identical(y,NA)){
        part <- factor(y['clusters'],levels=LETTERS[1:y['nbClusters']])
        plotLegend(as.numeric(table(part)/length(part)*100),adjustLegend=adjustLegend,parMean=parMean)
    }
## cat("\nD",screen())
    return(invisible())
}


##setMethod("plot",signature=c("ClusterLongData","missing"),.clusterLongData.plot)
ClusterLongData_plot <- function(x,y=NA,parTraj=parTRAJ(),parMean=parMEAN(),addLegend=TRUE,adjustLegend=-0.12,toPlot="both",
                     criterion=x["criterionActif"],nbCriterion=1000,xlab="Times",ylab=x["varNames"],closeScreenTraj=TRUE,...){

    if(any(is.na(y))){
        toPlot <- "traj"
        addLegend <- FALSE
    }else{}
    switch(EXPR=toPlot,
           "both"={
## cat("\nA",screen())
               listScreen <- split.screen(matrix(c(0,0.3,0.3,1,0,0,1,1),2))
               screen(2)
               ## ??? Liste des arguments a vérifier
               ClusterLongData_plotTrajMeans(x=x,y=y,parTraj=parTraj,parMean=parMean,addLegend=addLegend,adjustLegend=adjustLegend,xlab=xlab,ylab=ylab,...)

               screen(1)
               plotCriterion(as(x,"ListPartition"),criterion=criterion,nbCriterion=nbCriterion)
########################               if(closeScreenTraj){close.screen(1)}else{}

               if(closeScreenTraj){
                   close.screen(listScreen)
                   return(invisible())
               }else{
                   return(listScreen)
               }
           },
           "traj"={
## cat("\nB",screen())
               ClusterLongData_plotTrajMeans(x=x,y=y,parTraj=parTraj,parMean=parMean,addLegend=addLegend,adjustLegend=adjustLegend,xlab=xlab,ylab=ylab,...)
## cat("\nBB",screen())

           },
           "criterion"={
## cat("\nC",screen())
               plotCriterion(as(x,"ListPartition"),criterion=criterion,nbCriterion=nbCriterion)
## cat("\nCC",screen())
           }

    )
}

setMethod("plot",signature=c("ClusterLongData","numeric"),ClusterLongData_plot)
setMethod("plot",signature=c("ClusterLongData","Partition"),ClusterLongData_plot)


ClusterLongData_missing_plot <- function(x,y,parTraj=parTRAJ(),parMean=parMEAN(),toPlot="both",
                     criterion=x["criterionActif"],nbCriterion=1000,xlab="Times",ylab=x["varNames"],closeScreenTraj=TRUE,...){
     ClusterLongData_plot(x=x,y=NA,parTraj=parTraj,parMean=parMean,toPlot=toPlot, criterion=criterion,nbCriterion=nbCriterion,...)
}

setMethod("plot",signature=c("ClusterLongData","missing"),ClusterLongData_missing_plot)



ClusterLongData_plotTraj <- function(x,y,...){
   plot(x=x,y=y,parMean=parMEAN(type="n"),toPlot="traj",parTraj=parTRAJ(col="clusters"),...)
}
setMethod("plotTraj",signature=c("ClusterLongData","numeric"),ClusterLongData_plotTraj)


ClusterLongData_plotMeans <- function(x,y,...){
   plot(x=x,y=y,toPlot="traj",parTraj=parTRAJ(type="n"),...)
}
setMethod("plotMeans",signature=c("ClusterLongData","numeric"),ClusterLongData_plotMeans)



gald <- generateArtificialLongData <- function(
    nbEachClusters=50,time=0:10,varNames="V",
    meanTrajectories=list(function(t){0},function(t){t},function(t){10-t},function(t){-0.4*t^2+4*t}),
    personalVariation=function(t){rnorm(1,0,2)},
    residualVariation=function(t){rnorm(1,0,2)},
    decimal=2,percentOfMissing=0
){
    nbClusters <- length(meanTrajectories)
    if(length(nbEachClusters)==1){nbEachClusters <- rep(nbEachClusters,nbClusters)}else{}
    if(is.numeric(personalVariation)){eval(parse(text=paste("personalVariation <- function(t){rnorm(1,0,",personalVariation,")}",sep="")))}else{}
    if(length(personalVariation)==1){personalVariation <- rep(list(personalVariation),nbClusters)}else{}
    if(is.numeric(residualVariation)){eval(parse(text=paste("residualVariation <- function(t){rnorm(1,0,",residualVariation,")}",sep="")))}else{}
    if(length(residualVariation)==1){residualVariation <- rep(list(residualVariation),nbClusters)}else{}
    if(length(percentOfMissing)==1){percentOfMissing <- rep(percentOfMissing,nbClusters)}else{}
    nbTime <- length(time)
    idAll <- paste("i",1:(sum(nbEachClusters)),sep="")
    indivInCluster <- rep(1:nbClusters,times=nbEachClusters)

    traj <- matrix(NA,nrow=sum(nbEachClusters),ncol=nbTime)
    for (iIndiv in 1:nrow(traj)){
        traj[iIndiv,] <- meanTrajectories[[indivInCluster[iIndiv]]](time)+
                         personalVariation[[indivInCluster[iIndiv]]](time)+
                         apply(t(time),2,residualVariation[[indivInCluster[iIndiv]]])
    }
    traj <- round(traj,digits=decimal)


    for (iCluster in 1:nbClusters){
        nbVal <- nbTime*nbEachClusters[iCluster]
        while(sum(is.na(traj[indivInCluster==iCluster,]))/nbVal < percentOfMissing[iCluster]){
            randL <- floor(runif(1,cumsum(c(0,nbEachClusters))[iCluster]+1,cumsum(nbEachClusters)[iCluster]+1))
            randC <- floor(runif(1,1,nbTime+1))
            if(sum(!is.na(traj[randL,]))>1){traj[randL,randC]<-NA}else{}
        }
    }

    return(clusterLongData(traj,idAll=idAll,time=time,varNames=varNames))
}

cat("\n--------------------------------------------------------------------
------------------------- Class ClustLongData ----------------------
--------------------------------- Fin ------------------------------
--------------------------------------------------------------------\n")
