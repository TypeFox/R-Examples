meteo2STFDF <- function(obs,
                        stations,
                        obs.staid.time=c(1,2),
                        stations.staid.lon.lat=c(1,2,3),
                        crs=CRS(as.character(NA)),
                        delta=NULL
                        ) {
 
           
  ids <- unique(stations[,obs.staid.time[1]])
  
  time <- unique(obs[,obs.staid.time[2] ])
  time <- as.POSIXlt(sort(time))
  
  nt <- length(time)
  ns <- length(ids) # num os stations
  
  
  tempdf <- data.frame(rep(ids,nt),rep(time,each=ns) ) 
  names(tempdf) <- names(obs)[obs.staid.time]
  
  #   require(plyr)
  data <- join(tempdf,obs)
  
  data <- data[ order( data[,obs.staid.time[1] ]), ]
  data <- data[ order( data[,obs.staid.time[2] ]), ] # sort like 1st station 1st date, 2nd stations. 1st date ... check it carefuly 
  row.names(data) <- 1:length(data[,1])
  
  # system.time( merge(tempdf,obs, all=TRUE) )
  # system.time(join(tempdf,obs) )
  # join is 2 x faster
  ids <- data.frame(staid=ids)
  names(ids) <- names(stations) [ stations.staid.lon.lat[1] ]
  st <- join( ids, stations)
  names(st)[ stations.staid.lon.lat[2:3] ] <- c('lon', 'lat')
  coordinates(st) <-~ lon +lat
  st@proj4string <- crs
  
  data <- as.data.frame(data[,-obs.staid.time] )
  names(data)= names(obs)[-obs.staid.time]
  
  if (is.null(delta) && length(time)==1){
    endTime <- time + 86400
    stfdf <-STFDF(st, time , data, endTime)
  } else if (is.null(delta) && length(time)!=1) {
    stfdf <-STFDF(st, time , data)
  } else {
    endTime <- time + delta
    stfdf <-STFDF(st, time , data, endTime)
  }

  # count NAs per stations
  bools2 <- apply(matrix(stfdf@data[,1],
                         nrow=length(stfdf@sp),byrow=F), MARGIN=1,
                  FUN=function(x) sum(is.na(x)))
  # remove all NA
  stfdf1=stfdf[bools2!=nt,drop=F]
  
  row.names(stfdf@sp) <- 1:nrow(stfdf@sp)
  
  return(stfdf)
  
}
