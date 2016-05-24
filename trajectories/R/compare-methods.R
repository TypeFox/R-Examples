## compare tracks
setGeneric(
  name = "compare",
  def = function(tr1, tr2) standardGeneric("compare")
)

## get distances between 2 Track objects for each point in time where they overlap
## extend each track with these points
## create corresponding lines
## returns a difftrack object
compare.track <- function(tr1, tr2) {
  if (!requireNamespace("xts", quietly = TRUE))
    stop("package xts required for track comparison")
  if (!(xts::first(tr1@endTime) < xts::last(tr2@endTime) && xts::first(tr2@endTime) < xts::last(tr1@endTime)))
      stop("Time itervals don't overlap!")
  if (!identicalCRS(tr1, tr2))
      stop("CRS are not identical!")
  crs <- CRS(proj4string(tr1))
  track1.df <- cbind(as.data.frame(tr1)[c(coordnames(tr1), "time")])
  track2.df <- cbind(as.data.frame(tr2)[c(coordnames(tr2), "time")])  
  # intervals timestamps fall in
  ivs1 <- findInterval(track1.df$time, track2.df$time) 
  ivs2 <- findInterval(track2.df$time, track1.df$time)
  # find points and create new extended data frames
  newTrack1.df <- findPoints(track2.df, track1.df, ivs2)
  newTrack2.df <- findPoints(track1.df, track2.df, ivs1)
  # points on the original
  conns12 <- merge(newTrack2.df, track1.df, "time")
  conns21 <- merge(track2.df, newTrack1.df, "time")
  
  conns12 <- lineConnections(conns12, crs)
  conns21 <- lineConnections(conns21, crs)
    
  # extended tracks
  newTrack1 <- STIDF(SpatialPoints(cbind(newTrack1.df$x, newTrack1.df$y), crs), newTrack1.df$time, data.frame(1:nrow(newTrack1.df)))
  newTrack2 <- STIDF(SpatialPoints(cbind(newTrack2.df$x, newTrack2.df$y), crs), newTrack2.df$time, data.frame(1:nrow(newTrack2.df)))
  newTrack1 <- Track(newTrack1)
  newTrack2 <- Track(newTrack2)
  new("difftrack", track1 = newTrack1, track2 = newTrack2, conns1 = conns12, conns2 = conns21)
}

setMethod("compare", signature("Track"), compare.track)



## distances between 2 Tracks objects
setGeneric(
  name = "dists",
  def = function(tr1, tr2, ...) standardGeneric("dists")
)

## returns a matrix with the distance between each pair of tracks
dists.tracks <- function(tr1, tr2, f = mean, ...) {
  cols <- dim(tr1)[[1]] 
  rows <- dim(tr2)[[1]] 
  dists <- matrix(nrow=rows, ncol=cols) # matrix with NA's
  for (i in 1:cols) {
    for (j in 1:rows) {
      if (identical(f, frechetDist)) {
        dists[i,j] <- f(tr1[i], tr2[j])
      } else try({ ## try in case compare gives an error because tracks don't overlap in time
        difftrack <- compare(tr1[i], tr2[j])
        dists[i,j] <- f(c(difftrack@conns1@data$dists, difftrack@conns2@data$dists), ...)
      })
    }
  }
  dists
}
setMethod("dists", signature("Tracks"), dists.tracks)


## finds corresponding points for track1 on track2
findPoints <- function(tr1, tr2, ivs) {
  x <- tr2[,1]
  y <- tr2[,2]
  time <- tr2[,3]
  for (i in 1:nrow(tr1)) {
    if (!ivs[i] == 0 && !ivs[i] == nrow(tr2)) {
      iv <- ivs[i]
      tdiff1 <- tr1$time[i] - tr2$time[iv] # diff between timestamp and start of the interval it falls in
      tdiff2 <- tr2$time[iv+1] - tr2$time[iv] # diff between timestamps (calculated here because it often varies)
      ratio <- as.numeric(tdiff1)/as.numeric(tdiff2) 
      x1 <- tr2[iv,1] # segment coordinates
      y1 <- tr2[iv,2]
      x2 <- tr2[iv+1,1]
      y2 <- tr2[iv+1,2]
      x <- c(x, x1 + ratio * (x2 - x1)) #find point 
      y <- c(y, y1 + ratio * (y2 - y1))
      time <- c(time, tr1$time[i])
    }
  }
  newTrack <- data.frame(x, y, time)
  newTrack <- newTrack[!duplicated(newTrack),] # remove duplicates
  newTrack <- newTrack[order(newTrack$time),] # sort by timestamp
  newTrack
}

## creates SpatialLines 
lineConnections <- function(conns, crs) {
  Lines <- list()
  coords1 <- cbind(conns[,2], conns[,3])
  coords2 <- cbind(conns[,4], conns[,5])
  for (i in 1:nrow(conns)) {
    Lines <- c(Lines, list(Lines(Line(rbind(coords1[i,], coords2[i,])), ID = i)))
  }
  sl <- SpatialLines(Lines, crs)
  dists <- SpatialLinesLengths(sl)
  sl <- SpatialLinesDataFrame(sl, data.frame(time = conns$time, dists), FALSE)
  sl
}


## calculates the discrete frechet distance between two tracks
setGeneric(
  name = "frechetDist",
  def = function(track1, track2, ...) standardGeneric("frechetDist")
)

frechetDist.track <- function(track1, track2) {
  if (!requireNamespace("xts", quietly = TRUE))
    stop("package xts required for frechetDist comparison")
  if (!identicalCRS(track1, track2))
    stop("CRS are not identical!")
  dists <- spDists(track1@sp, track2@sp) #dists between all points
  dists[,1] <- cummax(dists[,1]) # cases where one of the trajectories is a point 
  dists[1,] <- cummax(dists[1,])
  for (i in 2:nrow(dists)) { # build rest of frechet distance matrix
    for (j in 2:ncol(dists)) {
      dists[i,j] <- max(dists[i,j], min(dists[i-1,j], dists[i-1,j-1], dists[i,j-1]))
    }
  }
  xts::last(xts::last(dists))
}

setMethod("frechetDist", signature("Track"), frechetDist.track)



## downsamples a track to the length of another one
setGeneric(
  name = "downsample",
  def = function(track1, track2, ...) standardGeneric("downsample")
)

# track1: track that will be downsampled
# track2: to the dimension of track2
downsample.track <- function(track1, track2) {
  if(!identicalCRS(track1, track2))
    stop("CRS are not identical!")
  if(dim(track1) == dim(track2))
    stop("Dimensions are euqal!")
  tr <- track1
  xy <- coordinates(track1)
  time <- index(track1@time)
  crs <- CRS (proj4string(track1))
  while(dim(track1) > dim(track2)) {
    d1 <- track1$distance # distances
    n <- length(d1) - 1 # number of segments between every second point
    xy1 <- cbind(head (xy, n), tail (xy, n))    
    d2.long <- head(d1, n) + tail(d1, n)
    xy.new <- list()
    for(i in 1:n) xy.new[[i]] <- rbind(head(xy, n)[i,], tail(xy, n)[i,])
    d2.short <- sapply (xy.new, function(x) spDists(x, longlat=TRUE)[1,2])
    remove <- which.min(d2.long - d2.short) + 1
    xy <- xy[- remove,]
    time <- time[- remove]
    stidf <- STIDF(SpatialPoints (xy, crs), time, data.frame(extraDat=rnorm(n)))
    tr  <- Track (stidf)
  }
  tr
}

setMethod("downsample", signature("Track"), downsample.track)
