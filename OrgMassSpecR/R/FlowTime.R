FlowTime <- function(internalDiameterMicrometers = 24, tubingLengthCentimeters = 45, 
                     flowRateMicrolitersPerMinute = 0.3) {
	internalDiameterMillimeters <- internalDiameterMicrometers / 1000
	tubingLengthMillimeters <- tubingLengthCentimeters * 10
	internalVolume <- pi * ((internalDiameterMillimeters / 2) ^ 2) * tubingLengthMillimeters   # mm^3 or uL
	time = (internalVolume / flowRateMicrolitersPerMinute) * 60   # sec
	message("time in seconds:")
	return(time)
}
