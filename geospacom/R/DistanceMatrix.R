DistanceMatrix <-
function(poly,id,unit=1000,longlat=TRUE,fun=distHaversine) {
  L <- length(poly)
  center <- list()
  center[L+1] = 1
  ids <- poly@data[[id]]
  if (!all.equal(unique(ids), ids)) {
    stop("The column for ids you chose is invalid('", id,
         "'), because not all its values are unique")
  }else{
    if(!length(ids)){
      stop("The specified id variable ('", id,"') is empty or not valid")
    }
  }
  ##  ids <- sort(ids)
  counter = 1
  for (i in ids) {
    center[counter] <- rgeos::gCentroid(poly[poly@data[[id]]==i,])
    counter <- (counter +1)
  }
  M <- matrix(0,L,L)
  colnames(M) <- ids
  rownames(M) <- ids
  if(longlat==TRUE){
    calibre <- 1/unit
    for (i in 1:(L-1)){
      for (j in (i+1):L){
        distance <- distm(center[[i]],center[[j]],fun)*calibre
        M[i, j] <- distance
        M[j, i] <- distance
      }
    }
  }else{
    calibre <- unit/10
    for (i in 1:(L-1)){
      for (j in (i+1):L){
        distance <- rgeos::gDistance(center[[i]],center[[j]])*calibre
        M[i, j] <- distance
        M[j, i] <- distance
      }
    }
  }  
  return(M)
}
