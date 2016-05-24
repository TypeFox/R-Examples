CalcBearing <-
function(initialLat,initialLong,finalLat,finalLong)
{
	initialLat<-initialLat/360*2*pi
	initialLong<-initialLong/360*2*pi
	finalLat<-finalLat/360*2*pi
	finalLong<-finalLong/360*2*pi
	atan2(sin(finalLong-initialLong)*cos(finalLat),cos(initialLat)*sin(finalLat)-sin(initialLat)*cos(finalLat)*cos(finalLong-initialLong))
}
