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

distFrechet <- function(Px,Py,Qx,Qy,timeScale=0.1,FrechetSumOrMax="max"){
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
    print(Mdist)
    print(Mfret)
    cat("\n\n Result =",Mfret[maxP,maxQ],"\n");
    return(Mfret[maxP,maxQ])
}




