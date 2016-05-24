DF2SpatialPointsDataFrame <- structure(function#change data.frame to SpatialPointsDataFrame
### This function modifies an object of class data.frame to one of class SpatialPointsDataFrame
(  x, ##<< data frame to be converted
   coords = c("x", "y"),##<< which columns are coordinates
   crs = sp::CRS("+init=epsg:28992") ##<< projection scheme
){
  #require(sp)
  sp::coordinates(x) <- ~x+y #c("x", "y")
  sp::proj4string(x) <- crs;
  return(x)
###  the new object of class SpatialPointsDataFrame
} , ex = function(){
  if (requireNamespace("sp", quietly = TRUE)) {
    data("meuse", package = "sp", envir = environment())
    meuseSP = DF2SpatialPointsDataFrame(meuse)
    
    sp::plot(meuseSP, asp = 1, cex = 4 * meuse$zinc/max(meuse$zinc),
         pch = 1, col = as.numeric(meuse$ffreq)+1 )
    data("meuse.riv", package = "sp", envir = environment())
    lines(meuse.riv)  
  } else {
    print("package sp must be installed for this example")
  }

  
})
