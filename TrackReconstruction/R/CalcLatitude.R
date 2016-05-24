CalcLatitude <-
function(initialLat, distance, bearing) #distance in meters
{
	initialLat<-initialLat/360*2*pi
	er<-6371000# earths radius in meters
	asin(sin(initialLat)*cos(distance/er)+cos(initialLat)*sin(distance/er)*cos(bearing)) #Latitude Rad
}
