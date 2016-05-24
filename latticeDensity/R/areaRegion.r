areaRegion <-
function(formLatticeOutput){
  #
  #  Computer the area of the bounding polygon
  #  and subtract the areas of the hole polygons
  #
require(splancs)
boundary <- formLatticeOutput$poly
hole.list <- formLatticeOutput$hole.list
number.holes <- length(hole.list)
hole.areas <- rep(NA,number.holes)
if(number.holes>0){
  for (i in 1:number.holes){
    hole.areas[i] <- areapl(hole.list[[i]])
  }
  areaRegion <- areapl(boundary) - sum(hole.areas)
}
if(number.holes==0){
  areaRegion <- areapl(boundary)
  }
return(areaRegion)
}

