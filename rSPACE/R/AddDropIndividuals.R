# Functions to add or remove individuals from a landscape

# Search for locations to place individuals -----------------
placeIndividuals <- function(map, use, N, buffer, wght=F){
   n <- length(use)

   if(wght==T){
     smp <- sample(n, n, prob=getValues(map)[use])
   } else {
     smp <- sample(n, n)
   }

   x <- coordinates(map)[use, 1]
   y <- coordinates(map)[use, 2]
   isLongLat <- as.numeric(grepl('+proj=longlat',proj4string(map)))

   use10 <- rep(0, n)
    use10[smp[1]] <- isLongLat

   USE = .C('sample_ind',
                  as.double(x),               #Longitude
                  as.double(y),               #Latitude
                  as.integer(N),              #Desired number of wolverines to place
                  as.double(buffer),          #Minimum distance between points
                  use = as.integer(use10),    #Indicator of whether to include or not
                  as.integer(n),              #Number of snow_points to check (maxid)
                  as.integer(smp - 1)         #Random order to use (shift for C indexing)
                )$use

   return((1:n)[USE==1])
}


# Add wolverines to simulation ----------------------------
addN<-function(dN, map, Parameters, wolv.df = NULL){
  habitatOK <- getValues(map) >= Parameters$habitat.cutoff
  nTypes <- length(Parameters$MFratio)

  if(is.null(wolv.df)){
   useLocations<-which(habitatOK)
   new.wolv<-lapply(1:nTypes,function(x){
                newInd<-placeIndividuals(map = map,
                          use    = useLocations,
                          N      = dN*Parameters$MFratio[x],
                          buffer = Parameters$buffer[x],
                          wght  = Parameters$wghts)
                return(useLocations[newInd])
              })
  } else {
   new.wolv<-lapply(1:nTypes, function(x){
              currentLocations <- wolv.df[wolv.df$type==x,]$locID
              currentLocations <- coordinates(map)[currentLocations, ]
              distances <- distanceFromPoints(map, currentLocations)
              territoryOK <- getValues(distances) > Parameters$buffer[x]

              useLocations <- which(habitatOK & territoryOK)
              newInd<-placeIndividuals(map = map,
                            use    = useLocations,
                            N      = round(dN*Parameters$MFratio[x]),
                            buffer = Parameters$buffer[x],
                            wght  = Parameters$wghts)

              useLocations[newInd]
              })
  }
  return(new.wolv)
}

# Remove wolverines from simulation ----------------------
dropN<-function(dN, map, Parameters, wolv.df){
  nTypes <- length(Parameters$MFratio)

  if(Parameters$wghts==T){
    drop.rows<-sample(nrow(wolv.df), dN, prob=getValues(map)[wolv.df$locID])
  } else {
    drop.rows<-sample(nrow(wolv.df),dN)
  }

  lost.wolv<-wolv.df[drop.rows, ]
  lost.wolv<-lapply(1:nTypes, function(x) lost.wolv$locID[lost.wolv$type==x])

  return(list(drop.rows, lost.wolv))
}

# Convert wolv list to data.frame -----------------------
wolv.dataframe<-function(wolv.list, map=NULL){
  nTypes <- length(wolv.list)
  nPerType <- sapply(wolv.list, length)

  wolv.df<-data.frame(type  = rep(1:nTypes, nPerType),
                      locID = unlist(wolv.list))
  if(!is.null(map)){
    xy<-coordinates(map)[wolv.df$locID, ]
    wolv.df$x<-xy[, 1]
    wolv.df$y<-xy[, 2]
    }

  return(wolv.df)
}
