HTSClusterUsersGuide <- function(view=TRUE)
#	Adapted from edgeR User's Guide
{
	f <- system.file("doc","HTSClusterUsersGuide.pdf",package="HTSCluster")
	if(view) {
		if(.Platform$OS.type == "windows") 
			shell.exec(f)
		else
			system(paste(Sys.getenv("R_PDFVIEWER"),f,"&"))
	}
	return(f)
}
