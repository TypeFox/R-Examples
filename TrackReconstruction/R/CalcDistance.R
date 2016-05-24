CalcDistance <-
function(initialLat, initialLong, finalLat, finalLong)
{
	initialLat<-initialLat/360*2*pi
	initialLong<-initialLong/360*2*pi
	finalLat<-finalLat/360*2*pi
	finalLong<-finalLong/360*2*pi
	acos(sin(initialLat)*sin(finalLat)+cos(initialLat)*cos(finalLat)*cos(finalLong-initialLong))*6371
}
