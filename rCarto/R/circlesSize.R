circlesSize <-
function(fdc,fixedNorm,pt,var,shareOfCircles,radiusMax,valueMax){
  # maximal extension of the map
  x1<-bbox(fdc)[1]
  y1<-bbox(fdc)[2]
  x2<-bbox(fdc)[3]
  y2<-bbox(fdc)[4]
  
  # choice of a circle normalization
  if (fixedNorm==FALSE){
    # normalization based on a share of the map occupied by circles
    # sum of the variable
    sc<-sum(as.numeric(pt[,var]),na.rm=TRUE)
    # maximum surface of the map
    sfdc<-(x2-x1)*(y2-y1)
    # circle radius(in map units)
    pt$varSize<-sqrt((pt[,var]*shareOfCircles*sfdc/sc)/pi)
  } else {
    # normalization based on a fixed maximum surface and value
    # maximum radius in cm
    radiusMax <- radiusMax/cm(1)
    # maximum radius in map units
    radiusMax <- xinch(radiusMax)
    # maximum surface in map units
    pgcs <- radiusMax*radiusMax*pi
    # surface in map units
    m <- pgcs/valueMax
    pt$varSize <- pt[,var]*m
    # circle radius (in map units)
    pt$varSize <- sqrt(pt$varSize/pi)
  }
  return (pt)
}
