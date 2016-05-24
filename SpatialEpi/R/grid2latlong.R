grid2latlong <-
function(input){

toradians <- atan(1)/45
radiusearth <- 0.5*(6378.2+6356.7)
sine51 <- sin( 51.5*toradians )


#-------------------------------------------------------------------------------
# If a Spatial Polygons
#-------------------------------------------------------------------------------
if(is(input)[1] == "SpatialPolygons"){
  for(i in 1:length(input@polygons)){
    # for all Polygons's in polygon
    for(j in 1:length(input@polygons[[i]]@Polygons)){
      # Convert coordinates
      new.coords <- as.matrix(cbind(
        input@polygons[[i]]@Polygons[[j]]@coords[,1]/(toradians*radiusearth*sine51),
        input@polygons[[i]]@Polygons[[j]]@coords[,2]/(toradians*radiusearth)
      ))
      colnames(new.coords) <- NULL
      rownames(new.coords) <- NULL
      
      # Update Polygons
      input@polygons[[i]]@Polygons[[j]]@coords <- new.coords
      input@polygons[[i]]@Polygons[[j]] <- Polygon(input@polygons[[i]]@Polygons[[j]])
    }
    # Update polygons
    input@polygons[[i]] <- Polygons(
      input@polygons[[i]]@Polygons,
      ID=input@polygons[[i]]@ID
    )	
  }
  output <- SpatialPolygons(input@polygons,proj4string=CRS("+proj=utm"))

  
#-------------------------------------------------------------------------------
# else return numeric
#-------------------------------------------------------------------------------
}else{
  output <- data.frame(cbind(
    x=input[, 1]/(toradians*radiusearth*sine51),
    y=input[, 2]/(toradians*radiusearth)
  ))		
}

return(output)	
}
