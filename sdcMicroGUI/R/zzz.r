.onAttach <- function(lib,pkg)
	{
  options("guiToolkit"="RGtk2")
  packageStartupMessage("\n--------\n")
  packageStartupMessage("you may start the graphical user interface by running 'sdcGUI()'\n")
}
