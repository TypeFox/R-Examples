###################################################################################
#' Start the SPOT GUI
#'
#' This function starts the graphical user interface for SPOT, which
#' is based on java. The .jar File started with this function is
#' in the directory where SPOT is installed. 
#'
#' Java runtime environment needs to be installed to use the GUI.
#'
#' @export
####################################################################################
spotGui<-function(){
	guiPath=find.package("SPOT");
	guiPath=file.path(guiPath, "GUI/spotGui.jar")
	system(paste("java -jar ",guiPath))
}
