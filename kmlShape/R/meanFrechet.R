############# Overview ##########
###
### [pour l'utilisateur] meanFrechet2 : calcule la moyenne de frechet entre 2 traj, avec les poids et le timeScale.
###    Elle elimine les manquantes, applique le timeScale puis appelle une fonction C. Au retour, elle corrige le timeScale
###
### [pour l'utilisateur] meanFrechet : calcule la moyenne de Frechet. Applique timeScale, vire les manquantes, puis appelle meanFrechetSample
###     ou meanFrechetHierarchic selon les cas.
###
### meanFrechetNbCouples : fusionne nbCouple trajectoire deux à deux, puis ajoute à la fin les trajectoires restantes.
###    -> Utilise un appel vers le C pour la fusion de deux traj
###    -> N'applique PAS timeScale et ne vire PAS les manquantes, car c'est fait par meanFrechet.
###
### meanFrechet2by2 : fusionne toutes les trajectoires jusqu'à en obtenir une seule.
###    Si n n'est pas une puissance de 2, fusionne des premiers couples en guise d'éliminatoire (avec meanFrechetNbCouples) pour obtenir une puissance de 2
###    Puis fusionne de proche en proche les couples deux à deux, jusqu'à ce qu'il ne reste qu'une trajectoire.
###    -> utilise meanFrechetNbCouples
###    -> N'applique PAS timeScale et ne vire PAS les manquantes, car c'est fait par meanFrechet.
###
### meanFrechetSample : méthode "all" (pour cela, sampleSize=NA) ou "sample". Ne gère pas les manquantes, ni le timeScale.
###    -> Défini l'ordre des trajectoires, puis appelles frechetMean2by2
###    -> N'applique PAS timeScale et ne vire PAS les manquantes, car c'est fait par meanFrechet.
###
### meanFrechetHierarchic : méthode hierarchical.
###    -> Utilise un appel vers le C pour la fusion de deux traj
###    -> N'applique PAS timeScale et ne vire PAS les manquantes, car c'est fait par meanFrechet.
###
### On ne s'occupe des manquantes et des timeScale que dans les deux premières, car se sont les seules appelables par l'utilisateur.
###

#meanFrechet2C <- function(Px,Py,Qx,Qy,weightPQ,FrechetSumOrMax){
#    lengthP <- length(Px)
#    lengthQ <- length(Qx)
#    way <- c("PQ","P","Q")
#    Mpath <- Mdist <- Mfret <- matrix(0,lengthP,lengthQ,dimnames=c(list(paste("P",1:lengthP,sep=""),paste("Q",1:lengthQ,sep=""))))
#
#    Mdist[1,1] <- dist(rbind(c(Px[1],Py[1]),c(Qx[1],Qy[1])))
#    Mfret[1,1] <- Mdist[1,1]
#    Mpath[1,1] <- "start"
#
#    for(i in 2:lengthP){
#        Mdist[i,1] <- dist(rbind(c(Px[i],Py[i]),c(Qx[1],Qy[1])))
#        Mfret[i,1] <- do.call(FrechetSumOrMax , list( Mfret[i-1,1] , Mdist[i,1] ) )
#        Mpath[i,1] <- "P"
#    }
#
#    for (j in 2:lengthQ){
#        Mdist[1,j] <- dist(rbind(c(Px[1],Py[1]),c(Qx[j],Qy[j])))
#        Mfret[1,j] <- do.call(FrechetSumOrMax , list( Mfret[1,j-1] , Mdist[1,j] ) )
#        Mpath[1,j] <- "Q"
#    }
#
 #   for(i in 2:lengthP){
  #      for (j in 2:lengthQ){
   #         Mdist[i,j] <- dist(rbind(c(Px[i],Py[i]),c(Qx[j],Qy[j])))
    #        movePQ <- Mfret[i-1,j-1]
     #       moveP <- Mfret[i-1,j]
#            moveQ <- Mfret[i,j-1]
#            Mfret[i,j] = do.call(FrechetSumOrMax , list( min(movePQ,moveP,moveQ) , Mdist[i,j] ) )
 #           Mpath[i,j] =  way[which.min(c(movePQ,moveP,moveQ))]
  #      }
   # }

#    print(Mdist)
 #   print(Mfret)
  #  print(Mpath)
#
#    i <- lengthP
 #   j <- lengthQ
  #  bestPath <- c(lengthP,lengthQ)
   # while(i>1||j>1){
    #    if(Mpath[i,j]=="Q"){
     #       j<-j-1
      #  }else{
       #     if(Mpath[i,j]=="P"){
        #        i<-i-1;
         #   }else{
          #      i<-i-1;
           #     j<-j-1;
            #}
#        }
 #       bestPath <- rbind(c(i,j),bestPath)
  #  }
   # colnames(bestPath)<-c("P","Q")
    #rownames(bestPath)<-NULL
   # print(bestPath)
#    return(
 #      data.frame(
  #        times=apply(cbind(Px[bestPath[,1]],Qx[bestPath[,2]]),1,function(x){weighted.mean(x,weightPQ)}),
   #       traj=apply(cbind(Py[bestPath[,1]],Qy[bestPath[,2]]),1,function(x){weighted.mean(x,weightPQ)}),
    #      weight=sum(weightPQ)
     #  )
 #   )
#}


meanFrechet2 <- function(Px,Py,Qx,Qy,timeScale=0.1,FrechetSumOrMax="sum",weightPQ=c(1,1)){
    missingsP <- is.na(Px)|is.na(Py); Px <- Px[!missingsP]*timeScale; Py <- Py[!missingsP]
    missingsQ <- is.na(Qx)|is.na(Qy); Qx <- Qx[!missingsQ]*timeScale; Qy <- Qy[!missingsQ]
    weightPQ <- weightPQ/sum(weightPQ)
    result <- .C("calcMeanFrechet", Px = as.numeric(Px), Py = as.numeric(Py), tailleP = as.integer(length(Px)), weightP=as.numeric(weightPQ[1]),
                     Qx = as.numeric(Qx), Qy = as.numeric(Qy),  tailleQ = as.integer(length(Qx)),  weightQ=as.numeric(weightPQ[2]),
 			   meanFrechetX = numeric(length(Px)+length(Qx)-2),meanFrechetY = numeric(length(Px)+length(Qx)-2),tailleMeanFrechet=as.integer(0),
                     sumOrMax=as.integer(FrechetSumOrMax=="sum"),PACKAGE="kmlShape")
#    result$times <- result$times/timeScale
    meanFrechet <- data.frame(times=(result$meanFrechetX/timeScale)[result$tailleMeanFrechet:1],traj=result$meanFrechetY[result$tailleMeanFrechet:1])
    return(meanFrechet)
}



###
### On considère que trajLong a nécessairement 4 colonnes qui sont, dans l'ordre :
###   - id
###   - times
###   - traj
###   - weight
### Si les noms des colonnes sont différents de id/times/traj/weight, on les change
### Pas de gestion des manquantes dans cette fonction, c'est fait avant.
### De même, trajLong[,1] DOIT être un integer allant de 1 à nbId

meanFrechetNbCouples <- function(trajLong,nbCouples,FrechetSumOrMax){
    names(trajLong) <- c("id","times","traj","weight")
    resultTraj <- trajLong[FALSE,]
    listId <- unique(trajLong[,1])
    maxId <-  max(listId)
    for(i in 1:nbCouples){
        Px <- trajLong[trajLong[,1]==listId[2*i-1],2]
        Py <- trajLong[trajLong[,1]==listId[2*i-1],3]
        Qx <- trajLong[trajLong[,1]==listId[2*i],2]
        Qy <- trajLong[trajLong[,1]==listId[2*i],3]
        weightPQ <- c(trajLong[trajLong[,1]==listId[2*i-1],4][1],trajLong[trajLong[,1]==listId[2*i],4][1])
        weightPQ <- weightPQ/sum(weightPQ)
        result <-  .C("calcMeanFrechet", Px = as.numeric(Px), Py = as.numeric(Py), tailleP = as.integer(length(Px)), weightP=as.numeric(weightPQ[1]),
                     Qx = as.numeric(Qx), Qy = as.numeric(Qy),  tailleQ = as.integer(length(Qx)),  weightQ=as.numeric(weightPQ[2]),
 			   meanFrechetX = numeric(length(Px)+length(Qx)-2),meanFrechetY = numeric(length(Px)+length(Qx)-2),tailleMeanFrechet=as.integer(0),
                     sumOrMax=as.integer(FrechetSumOrMax=="sum"),PACKAGE="kmlShape")

        meanFrechet <- data.frame(times=result$meanFrechetX[result$tailleMeanFrechet:1],traj=result$meanFrechetY[result$tailleMeanFrechet:1])
        resultTraj <- rbind(resultTraj,data.frame(id=2*i-1,meanFrechet,weight=sum(weightPQ)))
    }

    if(nbCouples*2<length(listId)){
        resultTraj <- rbind(resultTraj,trajLong[trajLong[,1] %in% listId[(nbCouples*2+1):length(listId)],])
    }else{}
    return(resultTraj)
}



meanFrechet2by2 <- function(trajLong,FrechetSumOrMax){
   listId <- unique(trajLong[,1])
   resultTraj <- trajLong

   #######
   ### Play-off : fusion des j=n-2^i couples
   nbOfPlayers <- length(listId)
   tournamentDeep <- floor(log(nbOfPlayers)/log(2))
   playoffs <- nbOfPlayers-2^tournamentDeep
   if(playoffs>0){
      resultTraj <- meanFrechetNbCouples(trajLong=trajLong,nbCouples=playoffs,FrechetSumOrMax=FrechetSumOrMax)
   }else{}

   for(iRound in tournamentDeep:1){
   # cat("### fusion\n")
      nbOfMatch <- 2^(iRound-1)
      resultTraj <- meanFrechetNbCouples(trajLong=resultTraj,nbCouples=nbOfMatch,FrechetSumOrMax=FrechetSumOrMax)
#    print(resultTraj)
   }

   # cat(" ### fini !\n")
   return(resultTraj[,2:3])
}



### Attention, ne gère pas les manquantes ??? A vérifier !
#meanFrechet2by2 <- function(trajLong,FrechetSumOrMax="max",timeScale=1){
#   trajLong[,2] <- trajLong[,2]*timeScale
#   result <- meanFrechet2by2C(trajLong=trajLong,FrechetSumOrMax=FrechetSumOrMax)
#   result$times <- result$times/timeScale
#   return(result)
#}



#meanFrechetAll <- function(listTraj,FrechetSumOrMax="max",distPoints=dist,timeScale=1,shuffle=TRUE){
#    if(shuffle){listTraj <- listTraj[sample(x=length(listTraj))]}else{}
#    return(meanFrechet2by2(listTraj=listTraj,FrechetSumOrMax=FrechetSumOrMax,distPoints=distPoints,timeScale=timeScale))
#}


meanFrechetSample <- function(trajLong,FrechetSumOrMax,shuffle=TRUE,sampleSize=NA){
    nbId <- length(unique(trajLong[,1]))
    if(is.na(sampleSize)){sampleSize <- nbId}else{}
    if(shuffle){
#       trajReordered <- data.frame()
 #      sampl <- sample(x=nbId,size=sampleSize)
  #     for(i in sampl){
   #         traj1 <- trajLong[trajLong[,1]==i,]
    #        trajReordered <- rbind(trajAll,traj1)
     #  }
        sampl <- sample(x=nbId,size=sampleSize)
        trajLong[,1] <- factor(trajLong[,1])
        trajLong[,1] <- factor(trajLong[,1], levels = levels(trajLong[,1])[sampl])
        trajReordered <- trajLong[order(trajLong[,1]),]
        trajReordered <- trajReordered[!is.na(trajReordered[,1]),]
        trajReordered[,1] <- as.integer(trajReordered[,1])
    }else{
       trajReordered <- trajLong[trajLong[,1] %in% unique(trajLong[,1])[1:sampleSize],]
    }
    return(meanFrechet2by2(trajLong=trajReordered,FrechetSumOrMax=FrechetSumOrMax))
}







meanFrechetHierarchic <- function(trajLong,FrechetSumOrMax,methodHclust="average"){
    nbId <- length(unique(trajLong[,1]))

    distH <- matrix(0,nbId,nbId)
    for(i in 1:(nbId-1)){for(j in i:nbId){
        distH[j,i] <- distFrechet(Px=trajLong[trajLong[,1]==i,2],Py=trajLong[trajLong[,1]==i,3],
             Qx=trajLong[trajLong[,1]==j,2],Qy=trajLong[trajLong[,1]==j,3], FrechetSumOrMax=FrechetSumOrMax)
    }}

    distH <- as.dist(distH)
    classifHierarchic <- hclust(as.dist(distH),methodHclust)

    listMoy <- list()
    orderMerge <- classifHierarchic$merge
    nbMerge <- nrow(orderMerge)
    for(i in 1:nbMerge){
        Px <- trajLong[trajLong[,1]==-orderMerge[i,1],2]
        Py <- trajLong[trajLong[,1]==-orderMerge[i,1],3]
        Qx <- trajLong[trajLong[,1]==-orderMerge[i,2],2]
        Qy <- trajLong[trajLong[,1]==-orderMerge[i,2],3]
        weightP <- trajLong[which.max(trajLong[,1]==-orderMerge[i,1]),4]
        weightQ <- trajLong[which.max(trajLong[,1]==-orderMerge[i,2]),4]

        result <- .C("calcMeanFrechet", Px = as.numeric(Px), Py = as.numeric(Py), tailleP = as.integer(length(Px)), weightP=as.numeric(weightP/(weightP+weightQ)),
                     Qx = as.numeric(Qx), Qy = as.numeric(Qy),  tailleQ = as.integer(length(Qx)),  weightQ=as.numeric(weightQ/(weightP+weightQ)),
 			   meanFrechetX = numeric(length(Px)+length(Qx)-2),meanFrechetY = numeric(length(Px)+length(Qx)-2),tailleMeanFrechet=as.integer(0),
                     sumOrMax=as.integer(FrechetSumOrMax=="sum"),PACKAGE="kmlShape")

        moy <- data.frame(-i,times=result$meanFrechetX[result$tailleMeanFrechet:1],traj=result$meanFrechetY[result$tailleMeanFrechet:1],weightPQ=weightP+weightQ)
        names(moy) <- names(trajLong)
        trajLong <- rbind(trajLong,moy)
    }
    return(data.frame(times=trajLong[trajLong[,1]==-nbMerge,2],traj=trajLong[trajLong[,1]==-nbMerge,3]))
}





meanFrechet <- function(trajLong,timeScale=0.1,FrechetSumOrMax="sum",aggregationMethod="all",shuffle=TRUE,sampleSize=NA,methodHclust="average"){
    if(ncol(trajLong)==3){
       trajLong$weight <- 1
    }else{
       if(ncol(trajLong)!=4){
           stop("[kmlShape:meanFrechet] The data.frame has to be (no choice) in the following format:
    - first column should be the individual indentifiant;
    - the second should be the times at which the measurement are made;
    - the third one should be the measurement;
    - the fourth (optional) can be the respective weight of each trajectories")
       }else{}
    }

    trajLong <- trajLong[!apply(trajLong,1,function(x)any(is.na(x))),]
    listId <- unique(trajLong[,1])
    nbId <- length(listId)

    if(nbId==1){
        result <- trajLong[,2:3]
    }else{
        trajLong[,1] <- as.integer(factor(trajLong[,1],labels=1:nbId))
        trajLong[,2] <- trajLong[,2]*timeScale

        if(aggregationMethod=="all"){
            result <- meanFrechetSample(trajLong=trajLong,FrechetSumOrMax=FrechetSumOrMax,shuffle=shuffle,sampleSize=NA)
        }else{
            if(aggregationMethod=="sample"){
                result <- meanFrechetSample(trajLong=trajLong,FrechetSumOrMax=FrechetSumOrMax,shuffle=shuffle,sampleSize=sampleSize)
            }else{
                if(aggregationMethod=="hierarchical"){
                    result <- meanFrechetHierarchic(trajLong=trajLong,FrechetSumOrMax=FrechetSumOrMax,methodHclust=methodHclust)
                }else{
                    stop("[meanFrechet] : this method is not implemented!")
                }
            }
        }
        result[,1] <- result[,1]/timeScale
    }
    names(result) <- c("times","traj")
    return(result)
}


### A été modifié, n'est plus opérationnel.
#meansFrechet <- function(trajLong,timeScale=0.1,aggregationMethod="all",sampleSize=0,FrechetSumOrMax="max",reroll=5,toPlot=FALSE){
#    indiv <- unique(trajLong[,1])
#    fDist <- function(i,curRes) {
#                distFrechet(Px = trajLong[trajLong[, 1] == i,2],
#                   Py = trajLong[trajLong[,1] == i, 3],
#                   Qx = curRes$times, Qy = curRes$traj,
#                   FrechetSumOrMax = FrechetSumOrMax, timeScale = timeScale)
#            }
#    bestDist <- +Inf
#    if(toPlot){matplotLong(trajLong,col="grey")}else{}
#    if(aggregationMethod=="hierarchical"){
#       for(i in 1:length(aggregationArgument)){
#           currentResult <- meanFrechet(trajLong=trajLong,timeScale=timeScale,aggregationMethod=aggregationMethod,shuffle=shuffle,
#                                        FrechetSumOrMax=FrechetSumOrMax,aggregationArgument=aggregationArgument[i])
#           if(toPlot){lines(currentResult,col=i+1,lwd=6)}else{}
#           distToMi <- sapply(indiv, fDist,currentResult)
#           currentDist <- mean(distToMi)
# #          print(criterion(Py=currentResult[,2],Qy=f(currentResult[,1]),times=currentResult[,1]))
#           if(currentDist<bestDist){
#               result <- currentResult
#               bestDist <- currentDist
#           }else{}
#       }

#    }else{
#       for(i in 1:reroll){
#           currentResult <- meanFrechet(trajLong=trajLong,timeScale=timeScale,aggregationMethod=aggregationMethod,shuffle=shuffle,
#                                        FrechetSumOrMax=FrechetSumOrMax,aggregationArgument=aggregationArgument)
#           if(toPlot){lines(currentResult,col=i+1,lwd=6)}else{}
#           distToMi <- sapply(indiv, fDist,currentResult)
#           currentDist <- mean(distToMi)
# #          print(criterion(Py=currentResult[,2],Qy=f(currentResult[,1]),times=currentResult[,1]))
#           if(currentDist<bestDist){
#               result <- currentResult
#               bestDist <- currentDist
#           }else{}
#       }
#    }

#    if(toPlot){
#        lines(result,col=1,lwd=6)
#        lines(result,col=2,lwd=3)
#    }else{}
#    return(result)
#}

