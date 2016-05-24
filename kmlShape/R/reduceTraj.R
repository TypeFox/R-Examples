 # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # reduceNbId  # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

reduceNbId <- function(id,trajWide,nbSenators=64,imputationMethod="linearInterpol"){
   if(any(is.na(trajWide))){trajWide <- imputation(trajWide,method=imputationMethod,lowerBound=NA,upperBound=NA)}else{}

   result <- kmeans(trajWide,centers=nbSenators,nstart=5)
   senatorsWeight <- as.integer(table(result$cluster))
   mySenators <- result$cluster
   senatorsWide <- result$centers

   reOrder <- sort(table(result$cluster),decreasing=TRUE)

   mySenators <- match(mySenators,names(reOrder))

   senatorsWeight <- as.integer(reOrder)
   names(senatorsWeight) <- paste("sen",1:nbSenators,sep="")
   senatorsWide <- senatorsWide[as.integer(names(reOrder)),]
   rownames(senatorsWide) <- names(senatorsWeight)
   mySenator <- data.frame(id=id,senators=paste("sen",mySenators,sep=""))

   return(list(mySenator=mySenator,senatorsWide=senatorsWide,senatorsWeight=senatorsWeight))
}





 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # reduceNbTimes # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

shortestDistanceToLines <- function(Mx,My,Ax,Ay,Bx,By){
    aire <- abs((By-Ay)*(Mx-Ax)-(Bx-Ax)*(My-Ay))
    return(  aire / sqrt((Bx-Ax)^2 + (By-Ay)^2))
}


findFarestPoint <- function(trajx,trajy){
    dmax <- 0
    index <- 1
    end <- length(trajx)

    if(end==2){
        index <- 1
        dmax <- 0
    }else{
        for(i in 2:(end-1)){
            d <- shortestDistanceToLines(Mx=trajx[i],My=trajy[i], Ax=trajx[1],Ay=trajy[1], Bx=trajx[end],By=trajy[end])
            if ( d > dmax ) {
                index <- i
                dmax <- d
            }else{}
        }
    }
    return(c(index=index,dmax=dmax))
}


DouglasPeuckerNbPoints <- function(trajx,trajy,nbPoints,spar=NA){
    missings <- is.na(trajx)|is.na(trajy)
    if(any(missings)){
        trajx <- trajx[!missings]
        trajy <- trajy[!missings]
    }else{}

    if(!is.na(spar)){trajy <- smooth.spline(trajx,trajy,spar=spar)[["y"]]}else{}

    result <- matrix(NA,nbPoints-1,4)
    colnames(result) <- c("firstPoint","lastPoint","middle","maxDist")
    result[1,] <- c(1,length(trajx),findFarestPoint(trajx,trajy))

    for(i in 2:(nbPoints-1)){
        if(max(result[,"maxDist"],na.rm=TRUE)==0){
            warning("[DouglasPeukerNbPoints] the simplified curve perfectly fit with the original with only ",i," points")
            break;
        }else{}
        lineToSplit <- which.max(result[,"maxDist"])
        firstPoint <- result[lineToSplit,1]
        lastPoint <- result[lineToSplit,2]
        middlePoint <- result[lineToSplit,3]

        result[lineToSplit,] <- c(firstPoint,middlePoint,findFarestPoint(trajx[firstPoint:middlePoint],trajy[firstPoint:middlePoint])+c(firstPoint-1,0))
        result[i,] <- c(middlePoint,lastPoint,findFarestPoint(trajx[middlePoint:lastPoint],trajy[middlePoint:lastPoint])+c(middlePoint-1,0))
    }
    x <- c(sort(result[,"firstPoint"]),length(trajx))
    return(data.frame(x=trajx[x],y=trajy[x]))
}



DouglasPeuckerEpsilon <- function(trajx,trajy,epsilon,spar=NA){
    missings <- is.na(trajx)|is.na(trajy)
    if(any(missings)){
        trajx <- trajx[!missings]
        trajy <- trajy[!missings]
    }else{}

    if(!is.na(spar)){trajy <- smooth.spline(trajx,trajy,spar=spar)[["y"]]}else{}

    farestPoint <- findFarestPoint(trajx,trajy)
    index <- farestPoint["index"]
    end <- length(trajx)
    if ( farestPoint["dmax"] > epsilon ) {
        recResults1 = DouglasPeuckerEpsilon(trajx[1:index],trajy[1:index], epsilon)
        recResults2 = DouglasPeuckerEpsilon(trajx[index:end],trajy[index:end], epsilon)

        resultTrajx = c(recResults1$x,recResults2$x[-1])
        resultTrajy = c(recResults1$y,recResults2$y[-1])
#        d = c(farestPoint["dmax"],recResults1$d,recResults2$d)
    } else {
        resultTrajx = c(trajx[1],trajx[end])
        resultTrajy = c(trajy[1],trajy[end])
#        d=numeric()
    }
    return(data.frame(x=resultTrajx,y=resultTrajy))
}


reduceNbTimes <- function(trajLong,nbPoints,spar=NA){
    if(ncol(trajLong)!=3){stop("[reduceNbTimes] The data.frame 'trajLong' has to be (no choice) in the following format:
    - first column should be the individual indentifiant;
    - the second should be the times at which the measurement are made;
    - the third one should be the measurement.")}else{}

    id <- unique(trajLong[,1])
    f <- function(i){
        traji <- trajLong[trajLong[,1]==i,]
        return(data.frame(i,DouglasPeuckerNbPoints(trajx=traji[,2],trajy=traji[,3],nbPoints=nbPoints,spar=spar)))
    }
    listDeTraj <- lapply(id,f)
    trajShort <- do.call(rbind,listDeTraj)
    colnames(trajShort) <- c("id","times","traj")
    return(trajShort)
}


 # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # reduceTraj  # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


reduceTraj <- function(myClds,nbSenators=NA,nbTimes=NA,spar=0.5,imputationMethod="linearInterpol"){
    if( is.na(nbSenators) & is.na(nbTimes) ){
       stop("[kmlShape:reduceTraj] Either 'nbSenators' or 'nbTimes' should be different from NA")
    }else{}

    nameObject <- deparse(substitute(myClds))
    ### Reduction of the number of individual
    if(!is.na(nbSenators)){
        myClds["reduceId"] <- TRUE

        if(!myClds["wideAvailable"]){convertTrajLongToWide(myClds,imputationMethod=imputationMethod)}else{}
#        show(myClds)


        resultSenat <- reduceNbId(id=myClds["id"],trajWide=myClds["trajWide"],nbSenators=nbSenators,imputationMethod=imputationMethod)
        myClds["senators"] <- reshapeWideToLong(data.frame(rownames(resultSenat$senatorsWide),resultSenat$senatorsWide),times=myClds["times"])
        myClds["mySenator"] <- resultSenat$mySenator
        myClds["senatorsWeight"] <- resultSenat$senatorsWeight
    }else{
        myClds["reduceId"] <- FALSE
        if(!myClds["longAvailable"]){convertTrajWideToLong(myClds)}else{}

        myClds["senators"] <- myClds["trajLong"]
        myClds["mySenator"] <- data.frame(myClds["id"],paste("sen",1:myClds["nbId"],sep=""))
        senatorsWeight <- as.integer(rep(1,length(myClds["id"])))
        names(senatorsWeight) <- paste("sen",1:myClds["nbId"],sep="")
        myClds["senatorsWeight"] <- senatorsWeight
    }

    ### Reduction of the number of time
    if(!is.na(nbTimes)){
        myClds["senators"] <- reduceNbTimes(myClds["senators"],nbPoints=nbTimes,spar=spar)
        myClds["reduceTimes"] <- TRUE
    }else{
        myClds["reduceTimes"] <- FALSE
    }

    myClds["senatorsAvailable"] <- TRUE

    assign(nameObject, myClds, envir = parent.frame())
    return(invisible())
}

