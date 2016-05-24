setAs("Move", "SpatialLinesDataFrame", function(from) {
  SpatialLinesDataFrame(SpatialLines(list(Lines(
    list(Line(from)),rownames(from@idData)
  )),proj4string = CRS(proj4string(from))), from@idData)
})
setAs("MoveStack", "SpatialLinesDataFrame", function(from) {
  SpatialLinesDataFrame(as(from, "SpatialLines"), from@idData)
})
setAs("MoveStack", "SpatialLines", function(from) {
  SpatialLines(
    mapply('Lines',
           lapply(lapply(
             lapply(split(from),as,"SpatialPoints")
             #      split(as(from,'SpatialPoints'), from@trackId)
             ,Line
           ), list)
           , rownames(from@idData)),proj4string = CRS(proj4string(from))
  )
})
