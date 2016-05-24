 # # # # # # # # # # # # # # # # # #
# # #     Class  definition     # # #
 # # #          clds           # # #
  # # # # # # # # # # # # # # # # #

### A mettre : champs de steps
Clds_validity <- function(object){
#    cat("**** validity clds <empty> ****\n")
    return(TRUE)
}

### A ajouter : nbId et nbCol
setClass(
    Class="Clds",
    representation=representation(
        steps="logical",
        id="factor",
        nbId="integer",
        nbCol="integer",
        trajWide="matrix",
        times="numeric",
        trajLong="data.frame",
        senators="data.frame",
        mySenator="data.frame",
        senatorsWeight="integer",
#        trajShort="data.frame",
        clustersSenators="factor",
        clusters="factor",
        trajMeans="data.frame"
   ),
    prototype=prototype(
        steps=logical(),
        id=factor(),
        nbId=integer(),
        nbCol=integer(),
        trajWide=matrix(,0,0),
        times=numeric(),
        trajLong=data.frame(),
        senators=data.frame(),
        mySenator=data.frame(),
        senatorsWeight=integer(),
 #       trajShort=data.frame(),
        clustersSenators=factor(),
        clusters=factor(),
        trajMeans=data.frame()
   ),
    validity=Clds_validity
)


cat("### Constructor ###\n")
cldsLong <- function(trajLong){
    id <- unique(trajLong[,1])
    return(new("Clds",
       steps=c(sourceWide=FALSE,wideAvailable=FALSE,longAvailable=TRUE,senatorsAvailable=FALSE,reduceId=FALSE,reduceTimes=FALSE,kmlShape=FALSE),
       id=as.factor(id),nbId=length(id),nbCol=integer(),
       trajLong=trajLong)
    )
}

### Deux formats pour trajWide :
###  - soit les id sont présent. Dans ce cas, trajWide est une matrice de trajectoire.
###  - soit les id sont absent. Dans ce cas, trajWide est un data.frame dont la première colonne est une matrice, les autres colonne les trajectoires

cat("### Constructor ###\n")
cldsWide <- function(trajWide,times,id){
    if(missing(id)){
       id <- factor(trajWide[,1])
       traj <- as.matrix(trajWide[,-1])
       if(missing(times)){times <- 1:(ncol(trajWide)-1)}else{}
    }else{
       id <- as.factor(id)
       traj <- trajWide
       if(missing(times)){times <- 1:ncol(trajWide)}else{}
    }
    rownames(traj) <- id
    return(new("Clds",
       steps=c(sourceWide=TRUE,wideAvailable=TRUE,longAvailable=FALSE,senatorsAvailable=FALSE,reduceId=FALSE,reduceTimes=FALSE,kmlShape=FALSE),
       id=id,
       trajWide=traj,
       times=times,
       nbId=nrow(traj),
       nbCol=ncol(traj)
    ))
}


cat("### Show ###\n")
Clds_show <- function(object){
    cat("   ~~~ Class:",class(object),"~~~ ")
    cat("\n ~ steps :\n")
    cat("     - Initial format of the data :",ifelse(object@steps["sourceWide"],"Wide","Long"),"\n")
    cat("     - Wide data available        :",object@steps["wideAvailable"],"\n")
    cat("     - Long data available        :",object@steps["longAvailable"],"\n")
    cat("     - Senators available         :",object@steps["senatorsAvailable"],"\n")
    cat("     - id reduction     :",object@steps["reduceId"],"\n")
    cat("     - times reduction  :",object@steps["reduceTimes"],"\n")
    cat("     - Use of kmlShape  :",object@steps["kmlShape"],"\n")
    if(object["wideAvailable"]){
       cat("\n ~ id : ")
       printLineShort(object@id)
       cat("\n ~ trajWide = [",object@nbId,"x",object@nbCol,"](limited to 5x10)\n",sep="")
       printMatrixShort(object@trajWide)
       cat("\n ~ times : ")
       catShort(object@times)
    }else{}
    if(object["longAvailable"]){
       cat("\n ~ trajLong = (limited to 10 lines)\n")
       printMatrixShort(object@trajLong,nRowToPrint=10)
    }else{}
    if(object["senatorsAvailable"]){
#reduceId"]|object["reduceTimes"]){
       cat("\n ~ senators (the five first):\n")
       printTrajLong(object@senators)
       cat("\n ~ mySenator :\n")
       printMatrixShort(object@mySenator,nRowToPrint=10)
       cat("\n ~ senatorsWeight : ")
       catShort(object@senatorsWeight)
    }else{}
    if(object["kmlShape"]){
       if(object["senatorsAvailable"]){
          cat("\n ~ clustersSenators : \n")
          printShort(object@clustersSenators)
       }else{}
       cat("\n ~ clusters : ")
       catShort(object@clusters)
       cat("\n ~ trajMeans : (the 25 first)\n")
       printTrajLong(object["trajMeans"],nRowToPrint=25)
    }
    cat("\n")
    return(invisible())
}

setMethod(f="show",signature="Clds",definition=Clds_show)


cat("### Getteur ###\n")
Clds_get <- function(x,i,j,drop){
    switch(EXPR=i,
       "steps"={return(x@steps)},
       "wideAvailable"={return(x@steps["wideAvailable"])},
       "longAvailable"={return(x@steps["longAvailable"])},
       "senatorsAvailable"={return(x@steps["senatorsAvailable"])},
       "reduceId"={return(x@steps["reduceId"])},
       "reduceTimes"={return(x@steps["reduceTimes"])},
       "kmlShape"={return(x@steps["kmlShape"])},
       "nbClusters"={return(max(x@trajMeans[,1]))},
       "id"={return(x@id)},
       "nbId"={return(x@nbId)},
       "nbCol"={return(x@nbCol)},
       "trajWide"={return(x@trajWide)},
       "times"={return(x@times)},
       "trajLong"={return(x@trajLong)},
       "senators"={return(x@senators)},
       "mySenator"={return(x@mySenator)},
       "senatorsWeight"={return(x@senatorsWeight)},
#       "trajShort"={return(x@trajShort)},
       "clustersSenators"={return(x@clustersSenators)},
       "clusters"={return(x@clusters)},
       "trajMeans"={return(x@trajMeans)},
       stop("[clds:get] ",i," is not a 'clds' slot")
    )
    return(invisible())
}
setMethod(f="[",signature="Clds",definition=Clds_get)

### Supprimer l'acces a id ? A trajW et trajL ?
cat("### Setteur ###\n")
Clds_set <- function(x,i,j,value){
    switch(EXPR=i,
       "steps"={x@steps<-value},
       "wideAvailable"={x@steps["wideAvailable"]<-value},
       "longAvailable"={x@steps["longAvailable"]<-value},
       "senatorsAvailable"={x@steps["senatorsAvailable"]<-value},
       "reduceId"={x@steps["reduceId"]<-value},
       "reduceTimes"={x@steps["reduceTimes"]<-value},
       "kmlShape"={x@steps["kmlShape"]<-value},
       "id"={x@id<-value},
       "nbId"={x@nbId<-value},
       "nbCol"={x@nbCol<-value},
       "trajWide"={x@trajWide<-value},
       "times"={x@times<-value},
       "trajLong"={x@trajLong<-value},
       "senators"={x@senators<-value},
       "mySenator"={x@mySenator<-value},
       "senatorsWeight"={x@senatorsWeight<-value},
 #      "trajShort"={x@trajShort<-value},
       "clustersSenators"={x@clustersSenators<-value},
       "clusters"={x@clusters<-value},
       "trajMeans"={x@trajMeans<-value},
       stop("[clds:set] ",i," is not a 'clds' slot")
    )
    validObject(x)
    return(x)
}
setMethod(f="[<-",signature="Clds",definition=Clds_set)


convertTrajLongToWide <- function(object,imputationMethod="linearInterpol"){
    nameObject <- deparse(substitute(object))
    if(object["wideAvailable"]){
        warning("[kmlShape:convertTrajLongToWide] trajWide already exists")
    }else{
        trajWide <- reshapeLongToWide(object["trajLong"])
        object["trajWide"] <- imputation(as.matrix(trajWide[,-1]),method=imputationMethod)
        object["times"] <- sort(unique(object["trajLong"][,2]))
        object["id"] <- as.factor(trajWide[,1])
        object["nbId"] <- nrow(trajWide)
        object["nbCol"] <- as.integer(ncol(trajWide)-1)
        object["wideAvailable"] <- TRUE
        assign(nameObject, object, envir = parent.frame())
    }
    return(invisible())
}



convertTrajWideToLong <- function(object){
    nameObject <- deparse(substitute(object))
    if(object["longAvailable"]){
        warning("[kmlShape:convertTrajWideToLong] trajLong already exists. It is recreated.")
    }else{}
    object["trajLong"] <- reshapeWideToLong(data.frame(object["id"],object["trajWide"]),times=object["times"])

    object["longAvailable"] <- TRUE
    assign(nameObject, object, envir = parent.frame())
    return(invisible())
}


