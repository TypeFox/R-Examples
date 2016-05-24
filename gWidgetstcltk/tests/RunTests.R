require(gWidgets)
require(tcltk)
options("guiToolkit"="tcltk")

if(as.numeric(.Tcl("info tclversion")) >= 8.5) {
  ## run tests only if we can
  files <- list.files(system.file("tests",package="gWidgets"),
                      pattern = "\\.R$",
                      full.names = TRUE)
  
  
  for(unitTest in files) {
    print(unitTest)
    source(unitTest)
  }
}
