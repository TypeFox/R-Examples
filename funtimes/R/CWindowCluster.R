CWindowCluster <- function(X, Alpha=NULL, Beta=NULL, Delta=NULL, Theta=0.8, p, w, s, Epsilon=1){
  T <- dim(X)[1]
  N <- dim(X)[2]
  nWindows <- length(seq(from=p*w, to=T, by=s*p))
  Clusters <- array(NA, dim=c(w, N, nWindows))
  ClustersWithinWindow <- array(NA, dim=c(nWindows, N))
  #Separate data into windows
  for (nw in 1:nWindows){
    #Apply CSlideCluster to each slide withing the window
    for (sl in 1:w){
      Clusters[sl, ,nw] <- CSlideCluster(X[c(((nw-1)*w*p+1+(sl-1)*p):((nw-1)*w*p+(sl)*p)),], Alpha=Alpha, Beta=Beta, Delta=Delta, Theta=Theta)
    }
    #Apply clusterization to CSlideCluster within window
    #count how many times each node was clustered together with each other node, frow w times
    E <- (sapply( 1:N, function(n) colSums(Clusters[,n,nw]==Clusters[,,nw])) >= Epsilon * w)
    TSclusters <- rep(NA, N)
    currentTS <- currentCluster <- 1
    TSclusters[currentTS] <- currentCluster
    Unclassified <- is.na(TSclusters)
    
    while(any(Unclassified)){
      Include <- CExpandWindowCluster(E[Unclassified,currentTS], E[Unclassified,Unclassified])
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
    ClustersWithinWindow[nw,] <- TSclusters
  }
  return(ClustersWithinWindow)
}
