cat("\n####################################################################
########################## Class LongData ##########################
############################# Creation #############################
####################################################################\n")

### Pas de trajectoire totalement vide => maxNA<length(time)

cat("### Definition ###\n")
LongData_validity <- function(object){
#    cat("**** validity LongData ****\n")
    if(length(object@idAll)==0&length(object@time)==0&length(object@varNames)==0&length(object@traj)==0){
    }else{
        if(any(c(length(object@idAll)==0,length(object@time)==0,length(object@varNames)==0,length(object@traj)==0))){
            stop("[LongData:validity]: at least one slot is empty")}else{}
        if(length(object@idFewNA)!=dim(object@traj)[1]){
            stop("[LongData:validity]: The number of id does not fit with the number of trajectories
  [LongData:validity]: length(idFewNA) =",length(object@idFewNA)," ; dim(traj)[1] =",dim(object@traj)[1])}else{}
        if(length(object@time)!=dim(object@traj)[2]){
            stop("[LongData:validity]: The number of time does not fit with the length of trajectories
  [LongData:validity]: length(time) =",length(object@time)," ; dim(traj)[2]=",dim(object@traj)[2])}else{}
        if(any(is.na(object@time))){
            stop("[LongData:validity]: There is some unknow times
  [LongData:validity]: is.na(time) =",is.na(object@time))}else{}
        if(!identical(object@time,sort(object@time))){
            stop("[LongData:validity]: time is not in increasing order
  [LongData:validity]: time =",object@time)}else{}
        if(any(duplicated(object@time))){
            stop("[LongData:validity]: Some time are duplicate
  [LongData:validity]: duplicated(time) =",duplicated(object@time))}else{}
        if(any(is.na(object@idAll))){
            stop("[LongData:validity]: Some idAll are NA
  [LongData:validity]: is.na(idAll) =",is.na(object@idAll))}else{}
        if(any(duplicated(object@idAll))){
            stop("[LongData:validity]: Some idAll are duplicate
  [LongData:validity]: duplicated(idAll) =",duplicated(object@idAll))}else{}
        if(any(dimnames(object@traj)[[1]]!=object@idFewNA,
               dimnames(object@traj)[[2]]!=paste("t",object@time,sep=""))){
            stop("[LongData:validity]: dimnames of traj is not correct
  [LongData:validity]: dimnames(traj) =",dimnames(object@traj),"
  [LongData:validity]: idFewNA =",object@idFewNA,"
  [LongData:validity]: paste('t',time) =",paste("t",object@time,sep=""))}else{}
        if(object@maxNA>=length(object@time)){
            stop("[LongData:validity]: maxNA is too high (trajectories with only NA are not trajectories)
  [LongData:validity]: maxNA =",object@maxNA," ; length(time) =",length(object@time))}else{}
    }
}

setClass(
    Class="LongData",
    representation=representation(
        idAll="character",
        idFewNA="character",
        time="numeric",
        varNames="character",
        traj="matrix",
        dimTraj="numeric",
        maxNA="numeric",
        reverse="matrix"
    ),
    prototype=prototype(
        idAll=character(),
        idFewNA=character(),
        time=numeric(),
        varNames=character(),
        traj=matrix(,0,0),
        dimTraj=numeric(),
        maxNA=numeric(),
        reverse=matrix(NA,2)
    ),
    validity=LongData_validity
)


cat("\n###################################################################
########################## Class LongData #########################
########################### Constructeur ##########################
###################################################################\n")


## buildObject <- function(traj,idAll,time,varNames,maxNA,reverse){
##     keepId <- apply(t(apply(traj,c(1,3),function(x){sum(is.na(x))}))<=maxNA,2,all)

##     ## Si on permet l'excusion globale, la formule est :
##     ## keepId <- apply(traj,1,function(x)(sum(is.na(x))<=maxNA))
##     traj <- traj[keepId,,,drop=FALSE]
##     idFewNA <- idAll
##     cat("VarNames=",idFewNA,"dimTraj",dim(traj),"\n")
##     dimnames(traj) <- list(idFewNA,paste("t",time,sep=""),varNames)
##     return(new("LongData",
##         idAll=as.character(idAll),
##         idFewNA=as.character(idFewNA),
##         time=time,
##         varNames=varNames,
##         traj=traj,
##         dimTraj=dim(traj),
##         maxNA=maxNA,
##         reverse=reverse)
##     )
## }



### Data.frame ou array en 2D
longData <- function(traj,idAll,time,timeInData,varNames,maxNA){
    if(missing(traj)){
       return(new("LongData"))
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
            stop("[LongData:constructor]: 'traj' should be either a data.frame, a matrix or an array")
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
    return(new("LongData",
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




cat("\n###################################################################
########################## Class LongData #########################
############################ Accesseurs ###########################
###################################################################\n")

cat("### Getteur ###\n")
LongData_get <- function(x,i,j,drop){
    switch(EXPR=i,
           "idAll"={return(x@idAll)},
           "idFewNA"={return(x@idFewNA)},
           "varNames"={return(x@varNames)},
           "time"={return(x@time)},
           "traj"={return(x@traj)},
           "dimTraj"={return(x@dimTraj)},
           "nbIdFewNA"={return(x@dimTraj[1])},
           "nbTime"={return(x@dimTraj[2])},
           "nbVar"={return(1)},
           "maxNA"={return(x@maxNA)},
           "reverse"={return(x@reverse)},
           stop("[LongData:get]:",i," is not a 'LongData' slot")
    )
}
setMethod("[","LongData",LongData_get)


### A priori, on n'a jamais besoin de modifier un LongData après sa création.
### ATTENTION : le set de ClusterLongData hérite directement de ListClustering
###    puisque set n'est pas défini pour LongData. Si on ajoute un set pour LongData,
###    il faut corriger le set de ClusterLongData

## cat("### Setteur ###\n")
## setMethod("[<-","LongData",
##     function(x,i,j,value){
##         traj <- x@traj
##         time <- x@time
##         varNames <- x@varNames
##         switch(EXPR=i,
##             "varNames"={varNames<-value},
##             "time"={time<-value},
##             "traj"={
##                 traj<-value
##                 if(length(dim(traj))==2){dim(traj) <- c(dim(traj),1)}else{}
##             },
##             stop("[LongData:setteur]: other slots can not be modified.")
##         )

##         dimnames(traj) <- list(x@idFewNA,paste("t",time,sep=""),varNames)
##         reverse <- x@reverse
##         dimnames(reverse) <- list(c("mean","sd"),varNames)
##         return(new("LongData",
##             idAll=x@idAll,
##             idFewNA=x@idFewNA,
##             time=time,
##             varNames=varNames,
##             traj=traj,
##             dimTraj=dim(traj),
##             maxNA=x@maxNA,
##             reverse=reverse)
##         )
##     }
## )



cat("\n###################################################################
########################## Class LongData #########################
############################# Affichage ###########################
###################################################################\n")

cat("### Method: 'show' pour LongData ###\n")
LongData_show <- function(object){
    cat("\n~ idAll       = [",length(object@idAll),"] ",sep="");catShort(object@idAll)
    cat("\n~ idFewNA     = [",object['nbIdFewNA'],"] ",sep="");catShort(object@idFewNA)
    cat("\n~ varNames    = [",object['nbVar'],"] ",sep="");catShort(object@varNames)
    cat("\n~ time        = [",object['nbTime'],"] ",sep="");catShort(object@time)
    cat("\n~ maxNA       = [",object['nbVar'],"] ",sep="");catShort(object@maxNA)
    cat("\n~ reverse     = [2x1]",sep="");
    cat("\n    - mean    =",object['reverse'][1,])
    cat("\n    - SD      =",object['reverse'][2,])
    cat("\n\n~ traj = [",object['nbIdFewNA'],"x",object['nbTime'],"] (limited to 5x10)  :\n",sep="")
    if(length(object@idFewNA)!=0){
        if(ncol(object@traj)>10){
            trajToShow <- as.data.frame(object@traj[,1:10])
            trajToShow$more <- "..."
        }else{
            trajToShow <- as.data.frame(object@traj)
        }
        if(nrow(object@traj)>5){
            print(trajToShow[1:5,])
            cat("... ...\n")
        }else{
            print(trajToShow)
        }
    }else{cat("   <no trajectories>\n")}
    return(invisible(object))
}

setMethod("show","LongData",
    definition=function(object){
        cat("\n   ~~~ Class: LongData ~~~")
        LongData_show(object)
    }
)


cat("### Method: 'print' pour LongData ###\n")
LongData_print <- function(x){
    object <- x
    cat("\n   ~~~ Class: LongData ~~~")
    cat("\n~ Class :",class(object))
    cat("\n\n~ traj = [",object['nbIdFewNA'],"x",object['nbTime'],"] (limited to 5x10)  :\n",sep="")
    print(object['traj'])
    cat("\n\n~ idAll = [",length(object@idAll),"]\n",sep="");print(object@idAll)
    cat("\n~ idFewNA = [",object['nbIdFewNA'],"]\n",sep="");print(object@idFewNA)
    cat("\n~ varNames = [",object['nbVar'],"]\n",sep="");print(object@varNames)
    cat("\n~ time = [",object['nbTime'],"]\n",sep="");print(object@time)
    cat("\n~ maxNA = [1]\n",sep="");print(object@maxNA)
    cat("\n~ reverse mean =\n");print(object['reverse'][1,])
    cat("\n~ reverse SD =\n");print(object['reverse'][2,])
    return(invisible(object))
}
setMethod("print","LongData",LongData_print)


setMethod("is.na", "LongData", function(x) FALSE) 



cat("\n###################################################################
########################## Class LongData #########################
############################## Various ############################
###################################################################\n")

LongData_scale <- function(x,center=TRUE,scale=TRUE){
    nameObject<-deparse(substitute(x))
    traj <- x@traj
    if(identical(center,TRUE)){center <- meanNA(traj)}else{}
    if(identical(scale,TRUE)){scale <- sdNA(as.numeric(traj))}else{}

    traj <- (traj-center)/scale
    x@reverse[1,] <- x@reverse[1,] + center*x@reverse[2,]
    x@reverse[2,] <- x@reverse[2,] * scale
    x@traj <- traj
    assign(nameObject,x,envir=parent.frame())
    return(invisible())
}

setMethod(f="scale",
    signature=c(x="LongData"),
    definition=LongData_scale
)


#.longData.scale2 <- function(x,center,scale){
#    traj <- x['traj']
#    if(missing(center)){center <- apply(traj,3,meanNA)}else{}
#    if(missing(scale)){scale <- apply(traj,3,sdNA)}else{}
#    traj <- sweep(traj,3,center)
#    traj <- sweep(traj,3,scale,FUN="/")
#    x@reverse[1,] <- center
#    x@reverse[2,] <- scale
#    x@traj <- traj
#    assign(nameObject,object,envir=parent.frame())
#    return(invisible())
#    x
#}
#setMethod(f="scale2",
#    signature=c(x="LongData"),
#    definition=.longData.scale2
#)


LongData_restoreRealData <- function(object){
    nameObject<-deparse(substitute(object))
    traj <- object@traj

    traj <- traj*object@reverse[2,1] + object@reverse[1,1]
    object@reverse[1,] <- 0
    object@reverse[2,] <- 1
    object@traj <- traj
    assign(nameObject,object,envir=parent.frame())
    return(invisible())
}
setMethod(f="restoreRealData",
    signature=c(object="LongData"),
    definition=LongData_restoreRealData
)


## gald <- generateArtificialLongData <- function(
##     nbEachClusters=50,time=0:10,varNames="V",
##     functionClusters=list(function(t){0},function(t){t},function(t){10-t},function(t){-0.4*t^2+4*t}),
##     constantPersonal=function(t){rnorm(1,0,2)},
##     functionNoise=function(t){rnorm(1,0,2)},
##     decimal=2,percentOfMissing=0
## ){
##     nbClusters <- length(functionClusters)
##     if(length(nbEachClusters)==1){nbEachClusters <- rep(nbEachClusters,nbClusters)}else{}
##     if(is.numeric(constantPersonal)){eval(parse(text=paste("constantPersonal <- function(t){rnorm(1,0,",constantPersonal,")}",sep="")))}else{}
##     if(length(constantPersonal)==1){constantPersonal <- rep(list(constantPersonal),nbClusters)}else{}
##     if(is.numeric(functionNoise)){eval(parse(text=paste("functionNoise <- function(t){rnorm(1,0,",functionNoise,")}",sep="")))}else{}
##     if(length(functionNoise)==1){functionNoise <- rep(list(functionNoise),nbClusters)}else{}
##     if(length(percentOfMissing)==1){percentOfMissing <- rep(percentOfMissing,nbClusters)}else{}
##     nbTime <- length(time)
##     idAll <- paste("i",1:(sum(nbEachClusters)),sep="")
##     indivInCluster <- rep(1:nbClusters,times=nbEachClusters)

##     traj <- matrix(NA,nrow=sum(nbEachClusters),ncol=nbTime)
##     for (iIndiv in 1:nrow(traj)){
##         traj[iIndiv,] <- functionClusters[[indivInCluster[iIndiv]]](time)+
##                          constantPersonal[[indivInCluster[iIndiv]]](time)+
##                          apply(t(time),2,functionNoise[[indivInCluster[iIndiv]]])
##     }
##     traj <- round(traj,digits=decimal)


##     for (iCluster in 1:nbClusters){
##         nbVal <- nbTime*nbEachClusters[iCluster]
##         while(sum(is.na(traj[indivInCluster==iCluster,]))/nbVal < percentOfMissing[iCluster]){
##             randL <- floor(runif(1,cumsum(c(0,nbEachClusters))[iCluster]+1,cumsum(nbEachClusters)[iCluster]+1))
##             randC <- floor(runif(1,1,nbTime+1))
##             if(sum(!is.na(traj[randL,]))>1){traj[randL,randC]<-NA}else{}
##         }
##     }

##     return(longData(traj,idAll=idAll,time=time,varNames=varNames))
## }


cat("\n-------------------------------------------------------------------
-------------------------- Class LongData -------------------------
------------------------------- Fin -------------------------------
-------------------------------------------------------------------\n")

