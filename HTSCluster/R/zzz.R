.onAttach <- function(libname, pkgname)
{
	if(.Platform$OS.type=="windows" && 
		.Platform$GUI=="Rgui") {

		winMenuAddItem("Vignettes", "HTSCluster",
		"shell.exec(system.file(\"doc\",\"HTSClusterUsersGuide.pdf\", package=\"HTSCluster\"))")

	}
}