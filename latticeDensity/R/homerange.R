homerange <- 
function(densityOut, percent = 0.95, output=FALSE){
#
      if(class(densityOut)!="densityOut"){
       stop("Should be the output from the function createDesity")}
  nodes <- densityOut$nodes
  poly <- densityOut$poly
  area <- densityOut$area
  z <- densityOut$probs
  cmstz <- cumsum(sort(z))
  count <- sum(cmstz <= (1-percent))
  ind <- (z>sort(z)[count])
  plot(nodes,cex=0.1)
  points(nodes[ind,],pch=19,cex=0.5)
  lines(rbind(poly,poly[1,]))
  #
  #  Compute proportion of total area in homerange
  #
  proportion <- sum(ind)/length(nodes[,1])
  prop.area <- proportion*area
  cat("Proportion of region in homerange = ", proportion, "\n")
  cat("Area of homerange = ", prop.area, "\n")
  if(output){
    return(cbind(nodes,ind))
    }
}

