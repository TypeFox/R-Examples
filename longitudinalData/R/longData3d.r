cat("\n####################################################################
######################### Class LongData3d #########################
############################# Creation #############################
####################################################################\n")

### Pas de trajectoire totalement vide => maxNA<length(time)

cat("### Definition ###\n")
LongData3d_validity <- function(object){
#    cat("**** validity LongData ****\n")
    if(length(object@idAll)==0&length(object@time)==0&length(object@varNames)==0&length(object@traj)==0){
    }else{
        if(any(c(length(object@idAll)==0,length(object@time)==0,length(object@varNames)==0,length(object@traj)==0))){
            stop("[LongData3d:validity]: at least one slot is empty")}else{}
        if(length(object@idFewNA)!=dim(object@traj)[1]){
            stop("[LongData3d:validity]: The number of id does not fit with the number of trajectories
  [LongData3d:validity]: length(idFewNA) =",length(object@idFewNA)," ; dim(traj)[1] =",dim(object@traj)[1])}else{}
        if(length(object@time)!=dim(object@traj)[2]){
            stop("[LongData3d:validity]: The number of time does not fit with the length of trajectories
  [LongData3d:validity]: length(time) =",length(object@time)," ; dim(traj)[2]=",dim(object@traj)[2])}else{}
        if(length(object@varNames)!=dim(object@traj)[3]){
            stop("[LongData3d:validity]: The number of variable does not fit with the width ot trajectories
  [LongData3d:validity]: length(varNames) =",length(object@varNames)," ; dim(traj)[3]=",dim(object@traj)[3])}else{}
        if(any(is.na(object@time))){
            stop("[LongData3d:validity]: There is some unknow times
  [LongData3d:validity]: is.na(time) =",is.na(object@time))}else{}
        if(!identical(object@time,sort(object@time))){
            stop("[LongData3d:validity]: time is not in increasing order
  [LongData3d:validity]: time =",object@time)}else{}
        if(any(duplicated(object@time))){
            stop("[LongData3d:validity]: Some time are duplicate
  [LongData3d:validity]: duplicated(time) =",duplicated(object@time))}else{}
        if(any(is.na(object@idAll))){
            stop("[LongData3d:validity]: Some idAll are NA
  [LongData3d:validity]: is.na(idAll) =",is.na(object@idAll))}else{}
        if(any(duplicated(object@idAll))){
            stop("[LongData3d:validity]: Some idAll are duplicate
  [LongData3d:validity]: duplicated(idAll) =",duplicated(object@idAll))}else{}
        if(any(dimnames(object@traj)[[1]]!=object@idFewNA,
               dimnames(object@traj)[[2]]!=paste("t",object@time,sep=""),
               dimnames(object@traj)[[3]]!=object@varNames)){
            stop("[LongData3d:validity]: dimnames of traj is not correct
  [LongData3d:validity]: dimnames(traj) =",dimnames(object@traj),"
  [LongData3d:validity]: idFewNA =",object@idFewNA,"
  [LongData3d:validity]: paste('t',time) =",paste("t",object@time,sep=""),"
  [LongData3d:validity]: varNames=",object@varNames)}else{}
        if(max(object@maxNA)>=length(object@time)){
            stop("[LongData3d:validity]: some maxNA are too high (trajectories with only NA are not trajectories)
  [LongData3d:validity]: maxNA =",object@maxNA," ; length(time) =",length(object@time))}else{}
    }
}

setClass(
    Class="LongData3d",
    representation=representation(
        idAll="character",
        idFewNA="character",
        time="numeric",
        varNames="character",
        traj="array",
        dimTraj="numeric",
        maxNA="numeric",
        reverse="matrix"
    ),
    prototype=prototype(
        idAll=character(),
        idFewNA=character(),
        time=numeric(),
        varNames=character(),
        traj=array(dim=c(0,0,0)),
        dimTraj=numeric(),
        maxNA=numeric(),
        reverse=matrix(NA,2)
    ),
    validity=LongData3d_validity
)


cat("\n###################################################################
########################## Class LongData #########################
########################### Constructeur ##########################
###################################################################\n")

### Data.frame ou array en 3D
longData3d <- function(traj,idAll,time,timeInData,varNames,maxNA){

    if(missing(traj)){
        return(new("LongData3d"))
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
            stop("[LongData3d:constructor]: 'traj' should be either a data.frame or an array")
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
    return(new("LongData3d",
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
LongData3d_get <- function(x,i,j,drop){
    switch(EXPR=i,
           "idAll"={return(x@idAll)},
           "idFewNA"={return(x@idFewNA)},
           "varNames"={return(x@varNames)},
           "time"={return(x@time)},
           "traj"={return(x@traj)},
           "dimTraj"={return(x@dimTraj)},
           "nbIdFewNA"={return(x@dimTraj[1])},
           "nbTime"={return(x@dimTraj[2])},
           "nbVar"={return(x@dimTraj[3])},
           "maxNA"={return(x@maxNA)},
           "reverse"={return(x@reverse)},
           stop("[LongData3d:get]:",i," is not a 'LongData' slot")
    )
}
setMethod("[","LongData3d",LongData3d_get)


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
##             stop("[LongData3d:setteur]: other slots can not be modified.")
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

cat("### Method: 'show' pour LongData3d ###\n")
LongData3d_show <- function(object){
    cat("\n~ idAll       = [",length(object@idAll),"] ",sep="");catShort(object@idAll)
    cat("\n~ idFewNA     = [",object['nbIdFewNA'],"] ",sep="");catShort(object@idFewNA)
    cat("\n~ varNames    = [",object['nbVar'],"] ",sep="");catShort(object@varNames)
    cat("\n~ time        = [",object['nbTime'],"] ",sep="");catShort(object@time)
    cat("\n~ maxNA       = [",object['nbVar'],"] ",sep="");catShort(object@maxNA)
    cat("\n~ reverse     = [2x",object['nbVar'],"]",sep="");
    cat("\n    - mean    =",object['reverse'][1,])
    cat("\n    - SD      =",object['reverse'][2,])
    cat("\n\n~ traj = [",object['nbIdFewNA'],"x",object['nbTime'],"x",object['nbVar'],"] (limited to 5x10x3)  :\n",sep="")
    if(length(object@idFewNA)!=0){
        for(iVar in 1:min(3,length(object@varNames))){
            cat("\n",object@varNames[iVar],":\n")
            if(ncol(object@traj)>10){
                trajToShow <- as.data.frame(object@traj[,1:10,iVar])
                trajToShow$more <- "..."
            }else{
                trajToShow <- as.data.frame(object@traj[,,iVar])
            }
            if(nrow(object@traj)>5){
                print(trajToShow[1:5,])
                cat("... ...\n")
            }else{
                print(trajToShow)
            }
        }
    }else{cat("   <no trajectories>\n")}
    return(invisible(object))
}

setMethod("show","LongData3d",
    definition=function(object){
        cat("\n   ~~~ Class: LongData3d ~~~")
        LongData3d_show(object)
    }
)


cat("### Method: 'print' pour LongData3d ###\n")
LongData3d_print <- function(x){
    object <- x
    cat("\n   ~~~ Class: LongData3d ~~~")
    cat("\n~ Class :",class(object))
    cat("\n\n~ traj = [",object['nbIdFewNA'],"x",object['nbTime'],"x",object['nbVar'],"] (limited to 5x10x3)  :\n",sep="")
    print(object['traj'])
    cat("\n\n~ idAll = [",length(object@idAll),"]\n",sep="");print(object@idAll)
    cat("\n~ idFewNA = [",object['nbIdFewNA'],"]\n",sep="");print(object@idFewNA)
    cat("\n~ varNames = [",object['nbVar'],"]\n",sep="");print(object@varNames)
    cat("\n~ time = [",object['nbTime'],"]\n",sep="");print(object@time)
    cat("\n~ maxNA = [",object['nbVar'],"]\n",sep="");print(object@maxNA)
    cat("\n~ reverse mean =\n");print(object['reverse'][1,])
    cat("\n~ reverse SD =\n");print(object['reverse'][2,])
    return(invisible(object))
}
setMethod("print","LongData3d",LongData3d_print)


setMethod("is.na", "LongData3d", function(x) FALSE)





cat("\n###################################################################
########################## Class LongData #########################
############################## Various ############################
###################################################################\n")

LongData3d_scale <- function(x,center=TRUE,scale=TRUE){
    nameObject<-deparse(substitute(x))
    traj <- x@traj
    if(identical(center,TRUE)){center <- apply(traj,3,meanNA)}else{}
    if(identical(scale,TRUE)){scale <- apply(traj,3,function(x){sdNA(as.numeric(x))})}else{}

    for (i in 1:x@dimTraj[3]){
        traj[,,i] <- (traj[,,i]-center[i])/scale[i]
    }
    x@reverse[1,] <- x@reverse[1,] + center*x@reverse[2,]
    x@reverse[2,] <- x@reverse[2,] * scale
    x@traj <- traj
    assign(nameObject,x,envir=parent.frame())
    return(invisible())
}

setMethod(f="scale",
    signature=c(x="LongData3d"),
    definition=LongData3d_scale
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


LongData3d_restoreRealData <- function(object){
    nameObject<-deparse(substitute(object))
    traj <- object@traj

    for (i in 1:object@dimTraj[3]){
        traj[,,i] <- traj[,,i]*object@reverse[2,i] + object@reverse[1,i]
    }
    object@reverse[1,] <- 0
    object@reverse[2,] <- 1
    object@traj <- traj
    assign(nameObject,object,envir=parent.frame())
    return(invisible())
}
setMethod(f="restoreRealData",
    signature=c(object="LongData3d"),
    definition=LongData3d_restoreRealData
)


## gald3d <- generateArtificialLongData3d <- function(
##     nbEachClusters=50,time=0:10,varNames=c("V","T"),
##     functionClusters=list(function(t){c(0,0)},function(t){c(10,10)},function(t){c(10-t,10-t)}),
##     constantPersonal=function(t){c(rnorm(1,0,2),rnorm(1,0,2))},
##     functionNoise=function(t){c(rnorm(1,0,2),rnorm(1,0,2))},
##     decimal=2,percentOfMissing=0#,clusterLongData=TRUE
## ){
##     nbClusters <- length(functionClusters)
##     if(length(nbEachClusters)==1){nbEachClusters <- rep(nbEachClusters,nbClusters)}else{}
##     if(is.numeric(constantPersonal)){eval(parse(text=paste("constantPersonal <- function(t){c(rnorm(1,0,",constantPersonal,"),rnorm(1,0,",constantPersonal,"))}",sep="")))}else{}
##     if(length(constantPersonal)==1){constantPersonal <- rep(list(constantPersonal),nbClusters)}else{}
##     if(is.numeric(functionNoise)){eval(parse(text=paste("functionNoise <- function(t){c(rnorm(1,0,",functionNoise,"),rnorm(1,0,",functionNoise,"))}",sep="")))}else{}
##     if(length(functionNoise)==1){functionNoise <- rep(list(functionNoise),nbClusters)}else{}
##     if(length(percentOfMissing)==1){percentOfMissing <- rep(percentOfMissing,nbClusters)}else{}
##     nbTime <- length(time)
##     nbVar <- length(varNames)
##     idAll <- paste("i",1:(sum(nbEachClusters)),sep="")
##     indivInCluster <- rep(1:nbClusters,times=nbEachClusters)

##     traj <- array(NA,dim=c(sum(nbEachClusters),nbTime,nbVar),dimnames=c(idAll,paste("t",time,sep=""),varNames))
##     for (iIndiv in 1:nrow(traj)){
## 	  traj[iIndiv,,] <- t(sapply(time,functionClusters[[indivInCluster[iIndiv]]])+constantPersonal[[indivInCluster[iIndiv]]](0)+sapply(time,functionNoise[[indivInCluster[iIndiv]]]))
##     }
##     traj <- round(traj,digits=decimal)

##     for (iCluster in 1:nbClusters){
##         nbVal <- nbTime*nbEachClusters[iCluster]
##         while(sum(is.na(traj[indivInCluster==iCluster,,]))/nbVal < percentOfMissing[iCluster]){
##             randL <- floor(runif(1,cumsum(c(0,nbEachClusters))[iCluster]+1,cumsum(nbEachClusters)[iCluster]+1))
##             randC <- floor(runif(1,1,nbTime+1))
##             randV <- floor(runif(1,1,nbVar+1))
##             if(sum(!is.na(traj[randL,,randV]))>1){traj[randL,randC,randV]<-NA}else{}
##         }
##     }
## #    if(clusterLongData){return(as.clusterLongData(traj,idAll=id,time=time,varNames=varNames))}else{
##     return(longData3d(traj,idAll=idAll,time=time,varNames=varNames))
## }

### allVarNames contient les nom des variables presentent dans un LongData.
### variable contient soit un nom de variable, soit le numero d'une variable.
### La fonction retourne le nom ET le numéro de la variable
varNumAndName <- function(variable,allVarNames){
    if(class(variable)=="character"){
        varName <- variable
        varNum <- c(1:length(allVarNames))[allVarNames %in% varName]
        if(length(varNum)==0){stop("[LongData3d:varNumAndName]: 'variable' is not a correct variable name
  [LongData3d:plod3d]: variable=",varName," is not in allVarNames=",allVarNames)}else{}
    }else{
        varNum <- variable
        varName <- allVarNames[varNum]
    }
    return(list(num=varNum,name=varName))
}


longDataFrom3d <- function(xLongData3d,variable){
    variable <- varNumAndName(variable,xLongData3d["varNames"])[[2]]
    selectVar <- xLongData3d["varNames"] %in% variable
    if(all(!selectVar)){stop("[LongData3d:longDataFrom3d] invalide variable names")}else{}
    idAll <- xLongData3d["idAll"]
    time <- xLongData3d["time"]
    traj <- xLongData3d["traj"][,,selectVar]
    traj <- rbind(traj,matrix(NA,nrow=length(idAll)-nrow(traj),ncol=ncol(traj),dimnames=list(idAll[!idAll %in% xLongData3d["idFewNA"]])))[idAll,]
    return(longData(traj=traj,
                    idAll=idAll,
                    time=time,
                    varNames=xLongData3d["varNames"][selectVar],
                    maxNA=xLongData3d["maxNA"][selectVar])
    )
}


longDataTo3d <- function(xLongData){
    idAll <- xLongData["idAll"]
    traj <- xLongData["traj"]
    time <- xLongData["time"]
    traj <- rbind(traj,matrix(NA,nrow=length(idAll)-nrow(traj),ncol=ncol(traj),dimnames=list(idAll[!idAll %in% xLongData["idFewNA"]])))[idAll,]
    dim(traj) <- c(dim(traj),1)
    return(longData3d(traj=traj,
                      idAll=idAll,
                      time=time,
                      varNames=xLongData["varNames"],
                      maxNA=xLongData["maxNA"])
    )
}





cat("\n-------------------------------------------------------------------
-------------------------- Class LongData -------------------------
------------------------------- Fin -------------------------------
-------------------------------------------------------------------\n")

