#####################################################
################### Distance Traj ###################
#####################################################

## distTraj <- function(x,y,method="euclidean",p=2){
##     power <- switch(method,
##                     "euclidean"=2,
##                     "manhattan"=1,
##                     "maximum"=Inf,
##                     "minkowski"=p,
##                     "canberra"=-1,
##                     "binary"=-2)
##     return(.C("distance",x=as.numeric(x),y=as.numeric(y),taille=as.integer(length(x)),
##             power=as.numeric(power),dist=as.numeric(0),NAOK=TRUE,PACKAGE="kmlShape")$dist)
## }



#####################################################
#################### Frechet max ####################
#####################################################

#####################
### Frechet max

distFrechetR <- function(Px,Py,Qx,Qy,timeScale=0.1,FrechetSumOrMax="sum"){
    missingsP <- is.na(Px)|is.na(Py);Px <- Px[!missingsP];Py <- Py[!missingsP]
    missingsQ <- is.na(Qx)|is.na(Qy);Qx <- Qx[!missingsQ];Qy <- Qy[!missingsQ]
    Px <- Px*timeScale
    Qx <- Qx*timeScale

    maxP <- length(Px)
    maxQ <- length(Qx)
    Mdist <- Mfret <- matrix(0,maxP,maxQ,dimnames=c(list(paste("P",1:maxP,sep=""),paste("Q",1:maxQ,sep=""))))
    for(i in 1:maxP){
        for (j in 1:maxQ){
            Mdist[i,j] <- dist(rbind(c(Px[i],Py[i]),c(Qx[j],Qy[j])))
            if(i == 1 && j == 1){Mfret[1,1] = Mdist[1,1]}
            if(i > 1 && j == 1){Mfret[i,1] = do.call(FrechetSumOrMax , list( Mfret[i-1,1] , Mdist[i,1] ) )}
            if(i == 1 && j > 1){Mfret[1,j] = do.call(FrechetSumOrMax , list( Mfret[1,j-1] , Mdist[1,j] ) )}
            if(i > 1  && j > 1){Mfret[i,j] = do.call(FrechetSumOrMax , list( min(Mfret[i-1,j],Mfret[i-1,j-1],Mfret[i,j-1]) , Mdist[i,j] ) )}
        }
    }
#    print(Mdist)
 #   print(Mfret)
  #  cat("\n\n Result =",Mfret[maxP,maxQ],"\n");
    return(Mfret[maxP,maxQ])
}


####################
### Frechet max selon l'article

distFrechetRec <- function(Px,Py,Qx,Qy,timeScale=0.1,FrechetSumOrMax="sum"){
    missingsP <- is.na(Px)|is.na(Py);Px <- Px[!missingsP];Py <- Py[!missingsP]
    missingsQ <- is.na(Qx)|is.na(Qy);Qx <- Qx[!missingsQ];Qy <- Qy[!missingsQ]
    Px <- Px*timeScale
    Qx <- Qx*timeScale

    maxP <- length(Px)
    maxQ <- length(Qx)
    Mdist <- matrix(0,maxP,maxQ,dimnames=c(list(paste("P",1:maxP,sep=""),paste("Q",1:maxQ,sep=""))))
    for(p in 1:maxP){
        for(q in 1:maxQ){
            Mdist[p,q] <- dist(rbind(c(Px[p],Py[p]),c(Qx[q],Qy[q])))
        }
    }
    ca <-matrix(-1,maxP,maxQ)
    rec <- function(i,j){
        if(ca[i,j] > -1){return(ca[i,j])}
        if(i == 1 && j == 1){ca[i,j] = Mdist[1,1]}
        if(i > 1 && j == 1){ca[i,1] = do.call( FrechetSumOrMax , list( rec(i-1,1) , Mdist[i,1] ) )}
        if(i == 1 && j > 1){ca[1,j] = do.call( FrechetSumOrMax , list( rec(1,j-1) , Mdist[1,j] ))}
        if(i > 1  && j > 1){ca[i,j] = do.call( FrechetSumOrMax , list( min(rec(i-1,j),rec(i-1,j-1),rec(i,j-1)) , Mdist[i,j]) )}
        return(ca[i,j]);
    }
    rec(maxP,maxQ)
}

########################
### Frechet en C
distFrechet <- function(Px,Py,Qx,Qy, timeScale=0.1,FrechetSumOrMax="sum"){
    missingsP <- is.na(Px)|is.na(Py);Px <- Px[!missingsP];Py <- Py[!missingsP]
    missingsQ <- is.na(Qx)|is.na(Qy);Qx <- Qx[!missingsQ];Qy <- Qy[!missingsQ]
    Px <- Px*timeScale
    Qx <- Qx*timeScale

        result <- .C("calcMatrixEuclidCumul",Px=as.numeric(Px),Py=as.numeric(Py),tailleP=as.integer(length(Px)),
                     Qx=as.numeric(Qx),Qy=as.numeric(Qy),tailleQ=as.integer(length(Qx)),
                     matDistEuclidCumul=numeric(length(Px)*length(Qx)),sumOrMax=as.integer(FrechetSumOrMax=="sum"),PACKAGE="kmlShape")
#    }else{
 #       if(FrechetSumOrMax=="sum"){
  #          result <- .C("distFrechetSum",Px=as.numeric(Px),Py=as.numeric(Py),Qx=as.numeric(Qx),Qy=as.numeric(Qy),
   #                      tailleP=as.integer(length(Px)),tailleQ=as.integer(length(Qx)),dist=as.numeric(0),NAOK=TRUE,PACKAGE="kmlShape")
    return(result$matDistEuclidCumul[length(Qx)*length(Px)])
}


##################################################
##################################################
##################################################


pathFrechet <- function(Px,Py,Qx,Qy,timeScale=0.1,FrechetSumOrMax="sum"){
    missingsP <- is.na(Px)|is.na(Py);Px <- Px[!missingsP];Py <- Py[!missingsP]
    missingsQ <- is.na(Qx)|is.na(Qy);Qx <- Qx[!missingsQ];Qy <- Qy[!missingsQ]
    Px <- Px*timeScale
    Qx <- Qx*timeScale

    result <- .C("calcPathFrechet",Px=as.numeric(Px), Py=as.numeric(Py),tailleP=as.integer(length(Px)),
        Qx=as.numeric(Qx),Qy=as.numeric(Qy),tailleQ=as.integer(length(Qx)),
        bestPathP=integer(length(Px)+length(Qx)-2),bestPathQ=integer(length(Px)+length(Qx)-2),tailleBestPath=as.integer(0),
        sumOrMax = as.integer(FrechetSumOrMax=="max"),PACKAGE="kmlShape")
#    }else{
#        if(FrechetSumOrMax=="sum"){
 #           result <- .C("pathFrechetSum",Px=as.numeric(Px),Py=as.numeric(Py),Qx=as.numeric(Qx),Qy=as.numeric(Qy),
  #                       tailleP=as.integer(length(Px)),tailleQ=as.integer(length(Qx)),dist=as.numeric(0),bestPath=integer((length(Px)+length(Qx))*2),pathLength=as.integer(0),NAOK=TRUE,PACKAGE="kmlShape")

    bestPath <- data.frame(P=result$bestPathP[result$tailleBestPath:1],Q=result$bestPathQ[result$tailleBestPath:1])+1
    return(bestPath)
}



#newQ <- numeric(maxP)
#for (p in 1:maxP){
#    newQ[p] <- mean(Q[fPath[,2][fPath[,1]==p]])
#}
#lines(newQ,col="green")
