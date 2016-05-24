cat("\n###################################################################
########################## Class Partition ########################
############################# Criterion ###########################
###################################################################\n")


matrix_qualityCriterion <- function(traj,clusters,imputationMethod="copyMean"){
    ## Si traj est un array, on colle les dimensions suivantes à la suite de la premiere

    if(nrow(traj)!=length(clusters)){
        stop("[qualityCriterion] the cluster and the number of trajectory should be the same.")
    }else{}

    clusters <- as.integer(clusters)
    nbIndiv <- nrow(traj)
    nbClusters <- max(clusters)
    if(max(nbClusters,na.rm=TRUE)==1){
        return(list(criters=c(
                Calinski.Harabatz=NA,Calinski.Harabatz2=NA,Calinski.Harabatz3=NA,Ray.Turi=NA,Davies.Bouldin=NA,
                BIC=NA,BIC2=NA,AIC=NA,AICc=NA,AICc2=NA,postProbaGlobal=NA,random=rnorm(1)),
                postProba=matrix(NA,nbIndiv,1),
                postProbaEachCluster=numeric()
           )
        )
    }else{}
    nbTime <- ncol(traj)

    if(any(is.na(traj))){trajImp <- imputation(traj,method=imputationMethod)}else{trajImp <- traj}
    trajImp <- matrix(as.numeric(trajImp),nrow=nrow(trajImp))       # Il arrive que traj soit une matrice d'entier, et ca coincerait...

    cls.attr <- cls.attrib(trajImp,clusters)

    #######################
    ### Calinski & Harabatz
    ###
    ### Selon Krzysztof     : C(k)=tB/tW*(n-1)/(n-k)
    ### Selon Milligan 1985 : C(k)=[tB/(k-1)]/[tW/(n-k)]=tB/tW*(n-k)/(k-1)

    varBetween <- bcls.matrix(cls.attr$cluster.center,cls.attr$cluster.size,cls.attr$mean)
    varWithin <- wcls.matrix(trajImp,clusters,cls.attr$cluster.center)
    traceBetween <- sum(diag(varBetween))
    traceWithin <- sum(diag(varWithin))
    calinski <- traceBetween/traceWithin*(nbIndiv-nbClusters)/(nbClusters-1)
    calinski2 <- traceBetween/traceWithin*(nbIndiv-1)/(nbIndiv-nbClusters)
    calinski3 <- traceBetween/traceWithin*(nbIndiv-nbClusters)/sqrt(nbClusters-1)
    if(is.na(calinski)){calinski<-NaN}

    #######################
    ### Ray & Turi
    ###
    ### Intra = moyenne des distances au carré entre un point et son centre
    ### Inter = plus petite distance au carré entre les centres.
    ### R(k)=Intra/Inter
    ### Un "grand" R dénote une mauvaise partition (grand Intra et/ou petit Inter)
    ###

    rayInter <- min(dist(cls.attr[[2]]))^2
    rayIntra <- mean((trajImp-cls.attr[[2]][clusters,])^2)
    ray <- as.numeric(rayIntra/rayInter)

    ## if(nrow(cls.attr[[2]])==partition@nbClusters){
    ##     rayInter <- +Inf
    ##     for(i in 1:partition@nbClusters){
    ##         distTrajI <- function(x){dist(rbind(x,cls.attr[[2]][i,]))}
    ##         rayInter <- min(apply(cls.attr[[2]][-i,,drop=FALSE],1,distTrajI),rayInter)
    ##     }
    ## }else{
    ##     rayInter <- NaN
    ## }

    ## rayIntra <- 0
    ## for (i in 1:nrow(traj)){
    ##   rayIntra <- rayIntra+dist(rbind(traj[i,],cls.attr[[2]][as.integer(clusters[i]),]))^2
    ## }
    ## rayIntra <- rayIntra/nrow(traj)
    ## ray <- as.numeric(rayInter/rayIntra)

    ##################
    ### Davies Bouldin
    ###
    ### Pour chaque cluster j, on définit MoyDistInt(j) une mesure des moyennes des distances internes de j (ou un diametre)
    ### Pour deux clusters j et j', DistExt(j,j') est une distance entre les clusters (exemple distance entre les centres de gravité)
    ### La 'proximite' entre j et j' vaut Proxi(j,j')=MoyDistInt(j)+MoyDistInt(j'))/DistExt(j,j')
    ### Si j et j' sont compacts et bien séparé, alors Proxi(j,j') sera petite.
    ###
    ### Ensuite, pour un cluster j, max(Proxi(j,j') donne sa pire proximité
    ### Au final, Davies and Bouldin est la moyenne des moins bonnes proximités de tous les clusters
    ### Un "grand" D denote une mauvaise partition (des proximités élevées)
    ###

    clsScat <- cls.scatt.data(trajImp,as.integer(clusters))
    davies <- as.numeric(clv.Davies.Bouldin(clsScat,"average","average"))

    ########################
    ### post probabilité, BIC, AIC
    ###
    ### On suppose la normalité des données. Cela permet de calculer la vraisemblance, puis tous les indices.
    ###

    preProba <- as.numeric(table(clusters))
    preProba <- preProba/sum(preProba)

    postProba <- vraisIndivXcluster <- matrix(,nbIndiv,nbClusters)

    ## Calcul de la vraisemblance
    moy <- matrix(,nbClusters,nbTime)
    for(i in 1:nbClusters){moy[i,] <- apply(traj[as.numeric(clusters)==i,,drop=FALSE],2,meanNA)}
    ecart <- sdcNA(as.numeric(traj-moy[clusters,]))
    for(i in 1:nbClusters){vraisIndivXcluster[,i] <- preProba[i]*apply(dnorm( t(traj),moy[i,],ecart ),2,prod,na.rm=TRUE)}
    postProba <- vraisIndivXcluster/apply(vraisIndivXcluster,1,sum)

    postProbaEachCluster <- rep(NA,nbClusters)
    for(i in 1:nbClusters){
        postProbaEachCluster[i] <- mean(postProba[as.numeric(clusters)==i,,drop=FALSE][,i])
    }
    postProbaGlobal <- mean(t(postProba)[as.numeric(clusters)+(0:(nbIndiv-1))*nbClusters,drop=FALSE])

    logVraisemblance <- sum(log(apply(vraisIndivXcluster,1,sum)))
    nbParam <- nbClusters*nbTime+1
    BIC  <- -2*logVraisemblance+nbParam*log(nbIndiv)
    BIC2 <- -2*logVraisemblance+nbParam*log(nbIndiv*nbTime)
    AIC  <- 2*nbParam-2*logVraisemblance
    AICc <- AIC+(2*nbParam*(nbParam+1))/(nbIndiv-nbParam-1)
    AICc2 <- AIC+(2*nbParam*(nbParam+1))/(nbIndiv*nbTime-nbParam-1)
#    entropie <- sum(apply(postProba,1,function(x) {ifelse(x==0,0,x*log(x))}))
 #   ICL <- BIC-2*entropie
  #  ICL2 <- BIC2-2*entropie

    return(list(criters=c(
                Calinski.Harabatz=calinski,Calinski.Harabatz2=calinski2,Calinski.Harabatz3=calinski3,Ray.Turi=-ray,Davies.Bouldin=-davies,
                BIC=-BIC,BIC2=-BIC2,AIC=-AIC,AICc=-AICc,AICc2=-AICc2,#entropie=entropie,ICL=-ICL,ICL2=-ICL2,
                postProbaGlobal=postProbaGlobal,random=rnorm(1)),
                postProba=postProba,
                postProbaEachCluster=postProbaEachCluster)
           )
}


setMethod("qualityCriterion",signature=c(traj="matrix",clusters="ANY",imputationMethod="ANY"),matrix_qualityCriterion)


array_qualityCriterion <- function(traj,clusters,imputationMethod="copyMean"){
    traj <- matrix(traj,nrow(traj))
    return(qualityCriterion(traj=traj,clusters=clusters,imputationMethod=imputationMethod))
}

setMethod("qualityCriterion",
          signature=c(traj="array",clusters="ANY",imputationMethod="ANY"),
          array_qualityCriterion
)



cat("####################################################################
########################## Class Partition #########################
############################# Creation #############################
####################################################################\n")

cat("### Definition ###\n")

Partition_validity <- function(object){
#    cat("**** validity Partition ****\n")
    if(!(length(object@nbClusters)==0&length(object@clusters)==0)){#not empty object
        if(any(c(length(object@nbClusters)==0,length(object@clusters)==0))){
            stop("[Partition:validity]: at least one slot is empty")}else{}
        if(object@nbClusters > MAX_CLUSTERS){
            stop("[Partition:validity]: More than ",MAX_CLUSTERS," clusters")}else{}
        if(!all(na.omit(object@clusters)%in%LETTERS[1:object@nbClusters])){
            stop("[Partition:validity]: Invalid clusters name, or clusters name out of range")}else{}
    }else{}
    return(TRUE)
}


### A priori, une partition peut contenir des clusters vides et n'est pas ordonnée.
### En pratique, kml et kml3d ne produisent que des ordonnées sans clusters vides.
### Le constructeur ne construit que des ordonnées sans vide, mais le vérificateur laisse tout de même la possibilité de créer 'manuellement' autre chose.
setClass(
   Class="Partition",
   representation=representation(
      nbClusters = "numeric",
      clusters = "factor",
      percentEachCluster="numeric",
      criterionValues="numeric",
      postProba="matrix",
      postProbaEachCluster="numeric",
      details="character" # algorithme, imputationMethod, convergenceTime, iteration, multiplicity
   ),
   prototype=prototype(
      nbClusters=numeric(),
      clusters=factor(),
      percentEachCluster=numeric(),
      criterionValues=numeric(),
      postProba=matrix(,0,0),
      postProbaEachCluster=numeric(),
      details=character()
   ),
   validity=Partition_validity
)

cat("\n####################################################################
########################## Class  Partition ########################
############################ Constructeur ##########################
####################################################################\n")


setMethod("partition",signature=c("missing","missing","missing"),function(){new("Partition")})

Partition_constructor <- function(clusters,traj,details=character()){
    ## Si clusters est numeric, il est transformé en LETTERS
    if(is.numeric(clusters)){
        if(max(clusters,na.rm=TRUE)>MAX_CLUSTERS){
            stop("[Partition:partition] the clusters should between 1 and ",MAX_CLUSTERS)
        }else{}
        clusters <- LETTERS[clusters]
    }else{}

    ## Vérification que clusters n'est QUE des LETTERS
    if(!all((clusters %in% LETTERS)|is.na(clusters))){
        stop("[Partition:partition] clusters should be either numeric or a vector of LETTERS")
    }else{}

    ## Ré ordonnancement + suppression des clusters vides
    tableClust <- table(clusters) ### Attention, cette table devient fausse apres le ré-ordonnancement des clusters
    nbClusters <- length(tableClust)
    nbIndiv <- length(clusters)

    clusters <- factor(clusters,
        levels=names(tableClust)[order(tableClust,-sapply(names(tableClust),function(x)which.max(clusters%in%x)),decreasing=TRUE)],
        labels=LETTERS[1:nbClusters]
    )

    ## Calcul des poucentages (= preProba)
    percentEachCluster <- as.numeric(table(clusters))
    percentEachCluster <- percentEachCluster/sum(percentEachCluster)

    ## Si un longData est fourni : calcul des post proba et des critères
    if(!missing(traj)){
        qualCriters <- qualityCriterion(traj,as.numeric(clusters))
    }else{
        qualCriters <- list(as.numeric(rep(NA,length(CRITERION_NAMES))),matrix(,0,0),as.numeric(NA))
    }
#    criters <- as.numeric(qualCriters[-c(1,2)])
 #   names(criters) <- CRITERION_NAMES

    return(new("Partition",nbClusters=nbClusters,clusters=clusters,percentEachCluster=percentEachCluster,
               criterionValues=qualCriters[[1]],details=details,postProba=qualCriters[[2]],postProbaEachCluster=qualCriters[[3]]))
}

setMethod("partition",signature=c(clusters="ANY",traj="missing",details="ANY"),Partition_constructor)

setMethod("partition",signature=c(clusters="ANY",traj="matrix",details="ANY"),Partition_constructor)
setMethod("partition",signature=c(clusters="ANY",traj="LongData",details="ANY"),
          function(clusters,traj,details){
              traj <- traj["traj"]
              return(partition(clusters,traj,details))
          }
)

Partition_constructor3d <- function(clusters,traj,details=character()){
    traj <- matrix(traj,nrow(traj))
    return(partition(clusters=clusters,traj=traj,details=details))
}

setMethod("partition",signature=c(clusters="ANY",traj="array",details="ANY"),Partition_constructor3d)
setMethod("partition",signature=c(clusters="ANY",traj="LongData3d",details="ANY"),
          function(clusters,traj,details){
              traj <- traj["traj"]
              return(partition(clusters,traj,details))
          }
)


cat("### Method : 'show' for partition ###\n") # Si on ajouter un titre a traj, on pourra afficher 'associate traj ='
Partition_show <- function(object){
    cat("   ~~~ Class : Partition ~~~ ")
    cat("\n ~ nbClusters           = ",object@nbClusters)
    cat("\n ~ percentEachCluster   = ",formatC(object@percentEachCluster,digits=2))
    cat("\n ~ qualities criterion:
   - Calinski & Harabatz                       =",object@criterionValues['Calinski.Harabatz'],"
   - Calinski & Harabatz modified by Kryszczuk =",object@criterionValues['Calinski.Harabatz2'],"
   - Calinski & Harabatz modified by Genolini  =",object@criterionValues['Calinski.Harabatz3'],"
   - Ray & Turie                        (opp.) =",object@criterionValues['Ray.Turi'],"
   - Davies & Bouldin                   (opp.) =",object@criterionValues['Davies.Bouldin'],"
   - BIC:  ln(L)-0.5xln(N)              (opp.) =",object@criterionValues['BIC'],"
   - BIC2: ln(L)-0.5xln(tN)             (opp.) =",object@criterionValues['BIC2'],"
   - AIC:  2ln(L)-2(2h)                 (opp.) =",object@criterionValues['AIC'],"
   - AICc: AIC + (2(2h)(2h+1))/(N-h-1)  (opp.) =",object@criterionValues['AICc'],"
   - AICc2: AIC + (2(2h)(2h+1))/(tN-h-1)(opp.) =",object@criterionValues['AICc2'],"
   - Overall post probilility                  =",object@criterionValues['postProbaGlobal'],"
   - random                                    =",object@criterionValues['random'])
    cat("\n ~ postProbaEachCluster = ",formatC(object@postProbaEachCluster,digits=2))
    cat("\n ~ convergenceTime      = ",as.integer(object@details["convergenceTime"]))
    cat("\n ~ multiplicity         = ",as.integer(object@details["multiplicity"]))
    cat("\n ~ imputationMethod     = ",object@details["imputationMethod"])
    cat("\n ~ algorithm            = ",object@details["algorithm"])
    cat("\n ~ clusters   : [",length(object@clusters),"]",sep="")
    if(length(object@nbClusters)!=0){
        for (iCluster in LETTERS[1:object@nbClusters]){
            toKeep <- iCluster==object@clusters
            cat("\n    ",iCluster," : [",sum(toKeep,na.rm=TRUE),"] ",sep="")
            catShort((1:length(object@clusters))[toKeep & !is.na(toKeep)])
        }
        cat("\n   <NA> : [",sum(is.na(object@clusters)),"] ",sep="")
        catShort((1:length(object@clusters))[is.na(object@clusters)])
        cat("\n")
    }else{
        cat("\n     <empty Partition>\n")
    }
    cat("\n ~ post probabilities (only the fist 5 lines):\n")
    if(nrow(object@postProba)>5){
        print(object@postProba[1:5,])
        cat("...")
    }else{
        print(object@postProba)
    }

    return(invisible(object))
}
setMethod(f="show",signature="Partition",definition=Partition_show)


cat("\n####################################################################
########################## Class Partition #########################
############################# Accesseurs ###########################
####################################################################\n")


cat("### Getteur ###\n")
setMethod("[","Partition",
    function(x,i,j,drop){
        switch(EXPR=i,
               "nbClusters"={return(x@nbClusters)},
               "clusters"={return(x@clusters)},
               "clustersAsInteger"={return(as.integer(x@clusters))},
               "percentEachCluster"={return(x@percentEachCluster)},
               "postProbaEachCluster"={return(x@postProbaEachCluster)},
               "postProba"={return(x@postProba)},
               "criterionValues"={return(x@criterionValues)},
               "convergenceTime"={return(as.integer(x@details["convergenceTime"]))},
               "multiplicity"={return(as.integer(x@details["multiplicity"]))},
               "imputationMethod"={return(as.integer(x@details["imputationMethod"]))},
               "algorithm"={return(x@details["algorithm"])},
               "details"={return(x@details)},
               if(i %in% CRITERION_NAMES){
                   return(x@criterionValues[i])
               }else{
                   if(i %in% names(x@details)){
                       return(x@details[i])
                   }else{
                       stop("[Partition:getteur]: there is not such a slot in Partition")
                   }
               }
        )
    }
)

cat("### Setteur ###\n")
setReplaceMethod("[","Partition",
    function(x,i,j,value){
        switch(EXPR=i,
               "multiplicity"={x@details["multiplicity"] <- as.character(value)},
               "convergenceTime"={x@details["convergenceTime"] <- value},
               if(i %in% c("clusters","nbClusters","percentEachCluster","algorithm","criterionNames","criterionValues")){
                   stop("[Partition:setteur]: ",i," is not entend to be change by the user.")
               }else{
                   stop("[Partition:setteur]: ",i," is not a 'Clustering' slot")
               }
        )
        validObject(x)
        return(x)
    }
)

setMethod("is.na", "Partition", function(x) FALSE) 


cat("\n####################################################################
########################## Class Partition #########################
############################### Autre ##############################
####################################################################\n")

## Comme le constructeur ordonne lui-même, cette fonction est sans doute obsolete

## .partition.ordered <- function(x){
##     clust <- x@clusters
##     tableClust <- table(clust)
##     nbClusters <- length(tableClust)
##     clusters <- factor(clust,
##         levels=names(tableClust)[order(tableClust,-sapply(names(tableClust),function(x)which.max(clust%in%x)),decreasing=TRUE)],
##         labels=LETTERS[1:nbClusters]
##     )

##     tab <- as.numeric(table(clusters))
##     x@percentEachCluster <- tab/sum(tab)
##     x@clusters <- clusters
##     x@nbClusters <- nbClusters
##     validObject(x)
##     return(x)
## }
## setMethod("ordered",signature="Partition",.partition.ordered)




LongData_qualityCriterion <- function(traj,clusters,imputationMethod="copyMean"){
    clust <- clusters["clustersAsInteger"]
    ##    resizePartition(traj,clusters)
    if(length(clust)!=traj["nbIdFewNA"]){
        clust <- clust[traj['idAll']%in%traj['idFewNA']]
    }else{}

    traj <- traj["traj"]
#    traj <- matrix(as.numeric(traj),nrow=nrow(traj))       # Il arrive que values soit une matrice d'entier, et ca coincerait...
    return(qualityCriterion(traj=traj,clusters=clust,imputationMethod=imputationMethod))
}

setMethod("qualityCriterion",
          signature=c(traj="LongData",clusters="Partition"),
          LongData_qualityCriterion
)


LongData3d_qualityCriterion <- function(traj,clusters,imputationMethod="copyMean"){
    clust <- clusters["clustersAsInteger"]
    ##    resizePartition(traj,clusters)
    if(length(clust)!=traj["nbIdFewNA"]){
        clust <- clust[traj['idAll']%in%traj['idFewNA']]
    }else{}

    traj <- traj["traj"]
#    traj <- matrix(as.numeric(traj),nrow=nrow(traj))       # Il arrive que values soit une matrice d'entier, et ca coincerait...
    return(qualityCriterion(traj=traj,clusters=clust,imputationMethod=imputationMethod))
}

setMethod("qualityCriterion",
          signature=c(traj="LongData3d",clusters="Partition"),
          LongData3d_qualityCriterion
)


initializePartition <- function(nbClusters,lengthPart,method="kmeans++",data){
    if(!missing(data) && !is.matrix(data)){
       data <- matrix(data,dim(data)[1])
    }else{}
    switch(method,
        "randomK"={
            part <- rep(NA,lengthPart)
            seeds <- sample(lengthPart,nbClusters)
            part[seeds] <- 1:nbClusters
        },
        "randomAll"={
            part <- floor(runif(lengthPart,1,nbClusters+1))       # Chaque individu recoit une affectation
            seeds <- sample(lengthPart,nbClusters)                # Puis on choisit k individus pour éviter les clusters vides.
            part[seeds] <- 1:nbClusters
        },
        "maxDist"={
            matrixDist <- as.matrix(dist(data))
            part <- rep(NA,lengthPart)
            seeds <- which(matrixDist==max(matrixDist,na.rm=TRUE),arr.ind=TRUE)[1,]
            part[seeds] <- 1:2
            nbSeeds <- 2
            while(nbSeeds<nbClusters){
                matrixDist[,seeds] <- 0
                nbSeeds <- nbSeeds+1
                seeds <- which.max(apply(matrixDist[!is.na(part),],2,min))[1]
                part[seeds] <- nbSeeds
            }
        },
        "kmeans++"={
            part <- rep(NA,lengthPart)
            seeds <- floor(runif(1,1,lengthPart+1))
            nbSeeds <- 1
            matrixDist <- matrix(NA,lengthPart,0)

            while(nbSeeds<nbClusters){
                newDist <- apply(data,1,function(x){dist(rbind(x,data[seeds[nbSeeds],]))})
                newDist[is.na(newDist)] <- max(newDist,na.rm=TRUE)
                matrixDist <- cbind(matrixDist,newDist)
                matrixDistMin <- apply(matrixDist,1,min)
                seeds <- c(seeds,sample(1:lengthPart,1,prob=matrixDistMin^2))
                nbSeeds <- nbSeeds+1
            }
            part[seeds] <- 1:nbClusters
        },
        ## kmeans++ en déterministe
        "kmeans+"={
            part <- rep(NA,lengthPart)
            seeds <- floor(runif(1,1,lengthPart+1))
            nbSeeds <- 1
            matrixDist <- matrix(NA,lengthPart,0)

            while(nbSeeds<nbClusters){
                matrixDist <- cbind(matrixDist,apply(data,1,function(x){dist(rbind(x,data[seeds[nbSeeds],]))}))
                matrixDistMin <- apply(matrixDist,1,min)
                seeds <- c(seeds,which.max(matrixDistMin))
                nbSeeds <- nbSeeds+1
            }
            part[seeds] <- 1:nbClusters
        },
        ## kmeans++ en supprimant le premier point choisi
        "kmeans--"={
            part <- rep(NA,lengthPart)
            seeds <- floor(runif(1,1,lengthPart+1))
            nbSeeds <- 1
            matrixDist <- apply(data,1,function(x){dist(rbind(x,data[seeds[nbSeeds],]))})
            seeds <- which.max(matrixDist)
            matrixDist <- matrix(NA,lengthPart,0)

            while(nbSeeds<nbClusters){
                newDist <- apply(data,1,function(x){dist(rbind(x,data[seeds[nbSeeds],]))})
		newDist[is.na(newDist)] <- max(newDist,na.rm=TRUE)
                matrixDist <- cbind(matrixDist,newDist)
                matrixDistMin <- apply(matrixDist,1,min)
                seeds <- c(seeds,sample(1:lengthPart,1,prob=matrixDistMin^2))
                nbSeeds <- nbSeeds+1
            }
            part[seeds] <- 1:nbClusters
        },
        ## kmeans++ en supprimant le premier point choisi
        "kmeans-"={
            part <- rep(NA,lengthPart)
            seeds <- floor(runif(1,1,lengthPart+1))
            nbSeeds <- 1
            matrixDist <- apply(data,1,function(x){dist(rbind(x,data[seeds[nbSeeds],]))})
            seeds <- which.max(matrixDist)
            matrixDist <- matrix(NA,lengthPart,0)

            while(nbSeeds<nbClusters){
                matrixDist <- cbind(matrixDist,apply(data,1,function(x){dist(rbind(x,data[seeds[nbSeeds],]))}))
                matrixDistMin <- apply(matrixDist,1,min)
                seeds <- c(seeds,which.max(matrixDistMin))
                nbSeeds <- nbSeeds+1
            }
            part[seeds] <- 1:nbClusters
        },
        stop("[PartitionInitialize] invalid initialization methods")
    )
    return(clusters=part)
}

#setMethod("initializePartition",signature=c("numeric","numeric","character","ANY"),initializePartition)

#initializePartitionArray <- function(nbClusters,lengthPart,method="kmeans++",data){
#    return(initializePartition(nbClusters=nbClusters,lengthPart=lengthPart,method=method,data=data))
#}

#setMethod("initializePartition",signature=c("numeric","numeric","character","array"),initializePartitionArray)



cat("\n--------------------------------------------------------------------
------------------------ Fin Class Partition -----------------------
--------------------------------------------------------------------\n")



