# Purpose        : shortest path between two points on WGS84 ellipsoid;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl);
# Contributions  : Dylan Beaudette (debeaudette@ucdavis.edu); 
# Status         : pre-alpha
# Note           : usefull for plotting airroutes or networks between geographically distant places;

geopath <- function(lon1, lon2, lat1, lat2, ID, n.points, print.geo = FALSE) {
    if(requireNamespace("fossil", quietly = TRUE)){
      ## Haversine Formula for Great Circle distance
      ## lon / lat = geographical coordinates on WGS84
      ## n.points = number of intermediate points
  
      p.1 <- matrix(c(lon1, lat1), ncol=2, dimnames=list(1,c("lon","lat")))  # source
      p.2 <- matrix(c(lon2, lat2), ncol=2, dimnames=list(1,c("lon","lat")))  # destination
      distc <- fossil::deg.dist(lat1=p.1[,2], long1=p.1[,1], lat2=p.2[,2], long2=p.2[,1])  # in km
      bearingc <- fossil::earth.bear(lat1=p.1[,2], long1=p.1[,1], lat2=p.2[,2], long2=p.2[,1])  # bearing in degrees from north
      # estimate the number of points based on the distance (the higher the distance, less points we need): 
      if(missing(ID)) { ID <- paste(ifelse(lon1<0, "W", "E"), abs(round(lon1,0)), ifelse(lat1<0, "S", "N"), abs(round(lat1,0)), ifelse(lon2<0, "W", "E"), abs(round(lon2,0)), ifelse(lat2<0, "S", "N"), abs(round(lat2,0)), sep="") }
      if(missing(n.points)) {
          n.points <- round(sqrt(distc)/sqrt(2), 0)
          }
      if(!is.nan(n.points)) { if(n.points>0) {
       pnts <- t(sapply(1:n.points/(n.points+1)*distc, FUN=fossil::new.lat.long, lat=p.1[,2], lon=p.1[,1], bearing=bearingc))[,c(2,1)] # intermediate points
      # some lines are crossing the whole globe (>180 or <-180 longitudes) and need to be split in two:
      if(is.matrix(pnts)){ if(max(LineLength(pnts, sum=FALSE))>100) {  
      breakp <- which.max(abs(pnts[,1]-c(pnts[-1,1], pnts[length(pnts[,1]),1])))
      pnts1 <- pnts[1:breakp,]
      pnts2 <- pnts[(breakp+1):length(pnts[,1]),]  
      routes <- Lines(list(Line(matrix(rbind(p.1, pnts1),ncol=2)), Line(matrix(rbind(pnts2, p.2),ncol=2))), ID=as.character(ID))
      } 
     
        else {
         routes <- Lines(list(Line(matrix(rbind(p.1, pnts, p.2),ncol=2))), ID=as.character(ID)) } 
  } } 
     # create SpatialLines:
     path <- SpatialLines(list(routes), CRS("+proj=longlat +datum=WGS84"))
     }
     if(print.geo==TRUE) {
     print(paste("Distance:", round(distc,1)))
     print(paste("Bearing:", round(bearingc,1)))
    } 
    return(path)
  }
}

# end of script;