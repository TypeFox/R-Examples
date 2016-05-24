################################################################################
# Converts SpatialLinesDataFrame to SpatialPointsDataFrame with points at
# the midpoints of the line segments (Jonathan Callahan).
SpatialLinesMidPoints <- function(sldf) {
  stopifnot(is.projected(sldf))
# define a function to get the midpoint of each Line object in each Lines
# object in slot(sldf, "lines"), copying out repeated roed of the data slot
  Lns <- slot(sldf, "lines")
  hash_lns <- sapply(Lns, function(x) length(slot(x, "Lines")))
  N <- sum(hash_lns)
  midpoints <- matrix(NA, ncol=2, nrow=N)
  Ind <- integer(length=N)
  ii <- 1
  for (i in 1:length(Lns)) {
    Lnsi <- slot(Lns[[i]], "Lines")
    for(j in 1:hash_lns[i]) {
      Ind[ii] <- i
      midpoints[ii,] <- getMidpoint(slot(Lnsi[[j]], "coords"))
      ii <- ii+1
    }
  }
  if (is(sldf, "SpatialLinesDataFrame")) {
    df0 <- slot(sldf, "data")[Ind,]
    df <- as.data.frame(cbind(df0, Ind))
  } else df <- data.frame(Ind=Ind)
# create a SpatialPointsDataFrame
  spdf <- SpatialPointsDataFrame(midpoints, data=df,
    proj4string=CRS(proj4string(sldf)))
  return(spdf)
}

getMidpoint <- function(coords) {
# calculate distances between points
  dist <- sqrt( (diff(coords[,1])^2 + (diff(coords[,2]))^2 ) )
# midpoint distance
  dist_mid <- sum(dist) / 2.0
# cumulative distances, starting with the first point
  dist_cum <- c(0,cumsum(dist))
# index of coordinates on either side of dist_mid
  end_index <- which(dist_cum > dist_mid)[1]
  start_index <- end_index-1
# calculate exact midpoint
  start <- coords[start_index,]
  end <- coords[end_index,]
  dist_remaining <- dist_mid - dist_cum[start_index]
  mid <- start + (end - start) * (dist_remaining / dist[start_index])
  return(mid)
}

