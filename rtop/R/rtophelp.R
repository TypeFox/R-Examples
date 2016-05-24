bbArea = function(bb) {
  xd = bb[[3]]-bb[[1]]
  yd = bb[[4]]-bb[[2]]
  abs(xd) * abs(yd)
}


summary.rtop = function(object,...) {
  summary.default(object,...)
}




adfunc = function(sampleVariogram, observations, amul) {
  if (is.null(sampleVariogram)) {
    if ("area" %in% names(observations)) {
      area = observations$area
    } else area = unlist(lapply(observations@polygons,FUN = function(poly) poly@area))
# alternative is variogram
  } else area = c(sampleVariogram$a1, sampleVariogram$a2)
  amax = max(area)
  amin = min(area)
  Rver = R.Version()
  if (as.numeric(Rver$major)*100 + as.numeric(Rver$minor) >= 214) {
    areas = axTicks(1,c(amin/5,amax*10,amul),usr = c(log10(amin/5)-1,log10(amax)+1), 
          log=TRUE, nintLog = Inf)
  } else {
    areas = axTicks(1,c(amin/5,amax*10,amul),usr = c(log10(amin/5)-1,log10(amax)+2), log=TRUE)
  }
  areas = areas[(min(which(areas > amin))-1):(max(which(areas < amax)) + 1)]
  areas
}




dfunc = function(sampleVariogram, observations, dmul) {
  if (is.null(sampleVariogram)) {
    dmax = sqrt(bbArea(bbox(observations)))/2
    dmin = min(dist(coordinates(observations)))
  } else if (is(sampleVariogram, "rtopVariogramCloud")) {
    dmax = max(sampleVariogram$dist) 
    dmin = min(sampleVariogram$dist)
  } else  {
    dmax = max(sampleVariogram$dist) 
    dmin = min(sampleVariogram$dist[sampleVariogram$np > 2])
  }
  Rver = R.Version()
  if (as.numeric(Rver$major)*100 + as.numeric(Rver$minor) >= 214) {
    dists = axTicks(1,c(dmin/5,dmax*10,dmul),usr = c(log10(dmin/5)-1,log10(dmax)+1), 
       log=TRUE, nintLog = Inf)
  } else {
    dists = axTicks(1,c(dmin/5,dmax*10,dmul),usr = c(log10(dmin/5)-1,log10(dmax)+2), log=TRUE)
  }
  dists[(min(which(dists > dmin))-1):(max(which(dists < dmax)) + 1)]
}

