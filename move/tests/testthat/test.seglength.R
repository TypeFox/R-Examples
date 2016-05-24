context("seglength")
test.seglength<-function(){
	a<-move(x=0:2,y=c(0,1,1),time=as.POSIXct(1:3, origin='1970-1-1'),proj=CRS('+proj=longlat +ellps=WGS84'))
#	require(geosphere)
	d<-distHaversine(coordinates(a)[-n.locs(a),],coordinates(a)[-1,] , 6378137)
	dd<-seglength(a)
	expect_equal(d,dd)
}


