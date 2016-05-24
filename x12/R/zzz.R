.onAttach <- function(...){
  #data(AirPassengers)
  packageStartupMessage("x12 is ready to use.")
  packageStartupMessage("Use the package x12GUI for a Graphical User Interface. \n")
  packageStartupMessage("It is advised to set the path to the X12 or X13 executables\n")
  packageStartupMessage("with x12path(validpath) or x13path(validpath)!\n")
  packageStartupMessage("---------------\n")
  packageStartupMessage("Suggestions and bug-reports can be submitted at: https://github.com/alexkowa/x12/issues")
}