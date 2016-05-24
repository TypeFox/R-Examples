CSlideCluster <- function(X, Alpha=NULL, Beta=NULL, Delta=NULL, Theta=0.8){
  #X is a matrix of time series, with one slide
  N <- dim(X)[2]
  if(is.null(Alpha) | is.null(Beta)){
    Alpha <- quantile(X)[2]-1.5*(quantile(X)[4]-quantile(X)[2])
    Beta <- quantile(X)[2]+1.5*(quantile(X)[4]-quantile(X)[2])
  }
  if(is.null(Delta)) {Delta <- 0.1*(Beta-Alpha)}
  TSclusters <- rep(NA, N)
  currentTS <- currentCluster <- 1
  TSclusters[currentTS] <- currentCluster
  Unclassified <- is.na(TSclusters)
  
  while(any(Unclassified)){
    Include <- CExpandSlideCluster(X[,currentTS], X[,Unclassified], Alpha, Beta, Delta, Theta)
    #Assign cluster number to the time series that were just groupped
    TSclusters[Unclassified][Include] <- currentCluster
    #Next cluster number will be:
    currentCluster <- currentCluster+1
    #Still unclassified time series
    Unclassified <- is.na(TSclusters)
    #Select the next time series to start the cluster with, and immediately
    #assign it to the new cluster
    currentTS <- which(Unclassified)[1]
    TSclusters[currentTS] <- currentCluster
    Unclassified <- is.na(TSclusters)
  }
  TSclusters
}
