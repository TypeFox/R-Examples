mutossGUI <- function(debugOutput=FALSE){
	#gui <- .jnew("mutoss.gui.MuTossGUI")
	.jcall("org/mutoss/gui/MuTossGUI", returnSig = "V", method="startGUI", debugOutput)
}

reportBug <- function(){
	.jcall("org/mutoss/gui/MuTossGUI", method="reportBug")
}
