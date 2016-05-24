mergeTrackObs<-function (sppdfInt, sppdfObs, obscol=NULL) {
  trc <- coordinates(sppdfInt)
  obs <- coordinates(sppdfObs)
  trckObs <- data.frame(trc, nObs = 0)
  for (i in 1:nrow(obs)) {
    distCal<-data.frame(matrix(rep(obs[i,],nrow(trc)),ncol=2,byrow=T),trc)
    distCal$dist<- sqrt((distCal[,1]-distCal[,3])^2+(distCal[,2]-distCal[,4])^2)
    idxtrc<-which(distCal$dist==min(distCal$dist))[1]
    if (is.null(obscol)){
      trckObs$nObs[idxtrc]<-trckObs$nObs[idxtrc]+1
    } else trckObs$nObs[idxtrc]<-trckObs$nObs[idxtrc]+sppdfObs@data[i,obscol]
  }
  trckObs <- data.frame(ID = 1:nrow(trckObs), trckObs)
  names(trckObs)[2:3] <- c("x", "y")
  coordinates(trckObs) <- ~x + y
  trckObs
}
