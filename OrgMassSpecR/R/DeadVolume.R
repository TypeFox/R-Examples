DeadVolume <- function(internalDiameterMicrometers = 24, tubeLengthCentimeters = 45) {
	internalDiameterMillimeters <- internalDiameterMicrometers / 1000
	tubeLengthMillimeters <- tubeLengthCentimeters * 10
	internalVolume <- pi * ((internalDiameterMillimeters / 2) ^ 2) * tubeLengthMillimeters   # mm^3 or uL
	message("internal volume in microliters:")
	return(internalVolume)
}
