degreeAxis = structure(function#axis with degrees
### add an axis with degree labels
(
  side, ##<< integer; see \link{axis}
  at, ##<< numeric; if missing, \link{axTicks} is called for nice values; see \link{axis}
  labels, ##<< character; if omitted labels are constructed with degree symbols, ending in N/S/E/W; in case of negative degrees, sign is reversed and S or W is added; see \link{axis}
  MyMap, ##<< optional map object to be passed
  ... ##<< optional arguments to \link{axis}
){

  degreeLabelsNS = function(x) {
	pos = sign(x) + 2
	dir = c("*S", "", "*N")
	paste(abs(x), "*degree", dir[pos])
  }
  degreeLabelsEW = function(x) {
	x <- ifelse(x > 180, x - 360, x)
	pos = sign(x) + 2
	if (any(x == -180))
		pos[x == -180] = 2
	if (any(x == 180))
		pos[x == 180] = 2
	dir = c("*W", "", "*E")
	paste(abs(x), "*degree", dir[pos])
  }
    USR = par('usr')
	
    if (missing(at)) {
        at = axTicks(side)
		#if (!missing(MyMap)) atLon = pretty(XY2LatLon(MyMap, X=USR[1:2],Y=USR[3])[,"lon"]) else atLon=at;
		#if (!missing(MyMap)) atLat = pretty(XY2LatLon(MyMap, X=USR[1],Y=USR[3:4])[,"lat"]) else atLat=at;
	}
	#browser();
    if (missing(labels)) {
        labels = FALSE
        if (side == 1 || side == 3) {
			if (!missing(MyMap)) atLon = XY2LatLon(MyMap, X=at,Y=USR[3])[,"lon"] else atLon=at;
            labels = parse(text = degreeLabelsEW(as.numeric(formatC(atLon, digits=5))))
        } else if (side == 2 || side == 4) {
		    if (!missing(MyMap)) atLat = XY2LatLon(MyMap, X=USR[1],Y=at)[,"lat"] else atLat=at;
            labels = parse(text = degreeLabelsNS(as.numeric(formatC(atLat, digits=5))))
		}
    }
  browser()
    ##note<< decimal degrees are used if variation is small, instead of minutes and seconds
    axis(side, at = at, labels = labels, ...)
### axis is plotted on current graph
}, ex = function(){
  xy = cbind(x = 2 * runif(100) - 1, y = 2 * runif(100) - 1)
  plot(xy,xlim=c(-1,1),ylim=c(-1,1))
  degreeAxis(1)
  degreeAxis(2, at = c(-1,-0.5,0,0.5,1))
})
