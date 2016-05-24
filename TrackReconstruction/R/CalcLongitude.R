CalcLongitude <-
function(initialLat, destinationLat, initialLong, distance, bearing) #distance in meters
{
	er<-6371000# earths radius in meters
	initialLat<-initialLat/360*2*pi
	initialLong<-initialLong/360*2*pi
	initialLong + atan2(sin(bearing)*sin(distance/er)*cos(initialLat),cos(distance/er)-sin(initialLat)*sin(destinationLat)) #Longitude Rad
}
