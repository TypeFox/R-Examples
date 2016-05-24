summary.spgeoOUT <- function(object, areanames = NA, ...) {
  suma <- paste(length(unique(object$identifier_in)), " species with ", dim(object$species_coordinates_in)[1], 
                " occurrence points and ", length(object$polygons), " input polygons.", sep = "")
  coords <- summary(object$species_coordinates)
  if (is.na(areanames) & length(object$areanam) != 0){
    areanames <- object$areanam
  }  
  if (is.na(areanames)){
    polys <- unlist(lapply(slot(object$polygons, "polygons"), function(x) slot(x, "ID")))
  }else{
    polys <- unique(as.character(data.frame(object$polygons)[, areanames]))
  }
  
  inf <- list(overall = suma, species_coordinates = coords, 
        polygons = polys, species_number_per_polygon = t(data.frame(Mean = mean(as.numeric(object$polygon_table)), 
            Median = median(as.numeric(object$polygon_table)), Max = max(object$polygon_table), Min = min(object$polygon_table))), 
        not_classified_samples = paste(dim(object$not_classified_samples)[1], " occurrences did not fall in any polygon", sep = ""))
    return(inf)
} 
