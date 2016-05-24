`qbbox` <-structure(function#computes bounding box
### The function qbbox computes a bounding box for the given lat,lon 
### points with a few additional options such as quantile boxes, additional margins, etc.
(
  lat, ##<<  latitude values
  lon, ##<<  longitude values
  TYPE = c("all", "quantile")[1], ##<< 
  margin = list(m=c(1,1,1,1), ##<< 
  TYPE = c("perc", "abs")[1]), ##<<  
  q.lat = c(0.1,0.9), ##<< 
  q.lon = c(0.1,0.9), ##<< 
  verbose=0 ##<<
){
 	if (TYPE == "all"){
 	  latR <- range(lat,na.rm=TRUE);
   	  lonR <- range(lon,na.rm=TRUE)
 	} else if (TYPE == "quantile"){
 	  latR <- quantile(lat, q.lat);
      lonR <- quantile(lon, q.lon);
 	}
 	if (!is.null(margin)){
 	  m <- margin$m;
 	  lat.center <- latR[1] + diff(latR)/2;
      lon.center <- lonR[1] + diff(lonR)/2;
      if (margin$TYPE == "perc"){
      	dlon <- c(-1,1)*(1+m[c(2,4)]/100)*diff(lonR)/2;
      	dlat <- c(-1,1)*(1+m[c(1,3)]/100)*diff(latR)/2;
      } else if (margin$TYPE == "abs"){
      	dlon <- c(-1,1)*(m[c(2,4)] + diff(lonR)/2);
      	dlat <- c(-1,1)*(m[c(1,3)] + diff(latR)/2);
      }
      lonR.margin <- lon.center + dlon;
      latR.margin <- lat.center + dlat;
      if (verbose>1) {
      	cat("old/new lon range:");print(lonR);print(lonR.margin);
      	cat("old/new lat range:");print(latR);print(latR.margin);
      }
      return(list(latR=latR.margin, lonR=lonR.margin))
 	}
 		
 	return(list(latR=latR, lonR=lonR))
### \item{latR }{latitude range}
###  \item{lonR }{longitude range}
 
 }, ex = function(){
  lat = 37.85 + rnorm(100, sd=0.001);
  lon = -120.47 + rnorm(100, sd=0.001);
  #add a few outliers:
  lat[1:5] <- lat[1:5] + rnorm(5, sd =.01);
  lon[1:5] <- lon[1:5] + rnorm(5, sd =.01);
  
  #range, discarding the upper and lower 10% of the data
  qbbox(lat, lon, TYPE = "quantile");
  #full range:
  qbbox(lat, lon, TYPE = "all");
  #add a 10% extra margin on all four sides:
  qbbox(lat, lon, margin = list(m = c(10, 10, 10, 10), TYPE = c("perc", "abs")[1]));
 
})


