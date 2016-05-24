# Purpose        : Derive distances to points;
# Maintainer     : Tomislav Hengl (tom.hengl@isric.org)
# Contributions  : ;
# Dev Status     : Pre-Alpha
# Note           : not recommended for big data;


setMethod("buffer.dist", signature(observations = "SpatialPointsDataFrame", predictionDomain = "SpatialPixelsDataFrame"), function(observations, predictionDomain, classes, width, ...){
  if(missing(width)){ width <- sqrt(areaSpatialGrid(predictionDomain)) }
  if(!length(classes)==length(observations)){ stop("Length of 'observations' and 'classes' does not match.") }
  s <- list(NULL)
  for(i in 1:length(levels(classes))){
    s[[i]] <- distance(rasterize(observations[which(classes==levels(classes)[i]),1]@coords, y=raster(predictionDomain)), width=width, ...)
  }
  s <- brick(s)
  s <- as(s, "SpatialPixelsDataFrame")
  s <- s[predictionDomain@grid.index,]
  return(s)
})