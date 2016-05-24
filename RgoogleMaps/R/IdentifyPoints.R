IdentifyPoints <-structure( function#identify points by clicking on map
### The user can try to identify lat/lon pairs on the map by clicking on them
( 
  MyMap, ##<< map object
  n=1,##<< the maximum number of points to locate.
  verbose = 0 ##<< level of verbosity
){
  cat("please identify ", n, "point(s) on the map with your mouse:\n")
  ret = locator(n)
  LatLon <- XY2LatLon(MyMap, ret$x,Y=ret$y)
  return(LatLon)
### the lat/lon coordinates of the chosen points are returned 
}, ex = function(){
  #The first step naturally will be to download a static map from the Google server. A simple example:
  
   #identifiy points:
   #IdentifyPoints(MyMap,5)
  
})
