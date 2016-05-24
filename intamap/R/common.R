coarsenGrid = function(object,coarse=2,offset = sample(c(0:(coarse-1)),2,replace=TRUE)) {
objCor = as(object,"SpatialPoints")
xss = unique(coordinates(objCor)[,1])
yss = unique(coordinates(objCor)[,2])
xi = c(1:(length(xss)/coarse))
xs = sort(xss)[(xi-1)*coarse+offset[1]+1]
yi = c(1:(length(yss)/coarse))
ys = sort(yss)[(yi-1)*coarse+offset[2]+1]
sel = coordinates(objCor)[,1] %in% xs & coordinates(objCor)[,2] %in% ys
newPoints = objCor[sel,]
if ("data" %in% names(getSlots(class(object)))) newPoints = SpatialPointsDataFrame(newPoints,data = object@data[sel,])
#  newPointDataFrame = as(object,"SpatialPointsDataFrame")[sel,]
gridded(newPoints) = TRUE
newPoints
}



 



commonArea = function(objecti,objectj) {
  bi = bbox(objecti)
  bj = bbox(objectj)
  iarea = bbArea(bi)
  jarea = bbArea(bj)
  sdim = sqrt((iarea+jarea)/2)
  bl = list()
  for (i in 1:2) bl[[i]] =  max(bi[[i]],bj[[i]])-sdim/1000
  for (i in 3:4) bl[[i]] =  min(bi[[i]],bj[[i]])+sdim/1000
  if (bl[[3]] >= bl[[1]] & bl[[4]] >= bl[[2]]) {
    larea = bbArea(bl)
    if (larea < 0.0001*max(iarea,jarea)) larea = 0.0001*max(iarea,jarea)
  } else {
    larea = 0
  }
#
  if (jarea< (0.0001*iarea)) jarea = 0.0001*iarea
  if (iarea< (0.0001*jarea)) iarea = 0.0001*jarea
  ilarea = larea/iarea
  jlarea = larea/jarea

#  print(bbox(objecti))
#  print(bbox(objectj))
#  print(bl)
#  cat(paste(iarea, jarea, larea, ilarea,jlarea,"\n"))

  return(list(ilarea,jlarea))
}


bbArea = function(bb) {
  xd = bb[[3]]-bb[[1]]
  yd = bb[[4]]-bb[[2]]
  abs(xd) * abs(yd)
}








SpatialDataFrame = function(spobj,data,...) {
  if (is(spobj,"SpatialPoints")) {
    SpatialPointsDataFrame(spobj,data = data,...)
  } else if (is(spobj,"SpatialGrid")) {
    SpatialGridDataFrame(spobj,data = data,...)
  } else if (is(spobj,"SpatialPixels")) {
    SpatialPixelsDataFrame(spobj,data = data,...)
  } else if (is(spobj,"SpatialLines")) {
    SpatialLinesDataFrame(spobj,data = data,...)
  } else if (is(spobj,"SpatialPolygons")) {
    SpatialPolygonsDataFrame(spobj,data = data,...)
  }
}
