summary.spgeoIN <- function(object, areanames = NA, ...) {
    
  sum <- paste(length(unique(object$identifier)), " species with ", dim(object$species_coordinates)[1], " occurrence points and ", 
               length(object$polygons), " input polygons.", sep = "")
  
  coords <- summary(object$species_coordinates)
  if (is.na(areanames) & length(object$areanam) != 0){
    areanames <- object$areanam
  }    
  if (is.na(areanames)){
    polys <- unlist(lapply(slot(object$polygons, "polygons"), function(x) slot(x, "ID")))
  }else{
    polys <- unique(data.frame(object$polygons)[, areanames])
  }
  
  suma <- list(overall = sum, species_coordinates = coords, polygons = polys)
  return(suma)
} 
