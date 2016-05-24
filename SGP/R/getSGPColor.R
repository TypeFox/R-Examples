`getSGPColor` <-
function(
	sgp,
	type) {

	if (toupper(type)=="PRESENTATION") {
		return(colorRampPalette(c("red","yellow3","springgreen","royalblue"))(99)[sgp])
	}
	if (toupper(type)=="PRINT") {
		return(colorRampPalette(c("red","yellow3","springgreen","royalblue"))(99)[sgp])
	}
} ### END getSGPColor function
