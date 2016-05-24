.onAttach <- function(lib, pkg){
  info <- packageDescription("FIACH")
  packageStartupMessage(
    paste('FIACH ',info$Version,":"," For the Graphical User Interface enter GUI() into the console", sep="")
 )
}


