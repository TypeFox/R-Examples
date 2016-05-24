#' @importFrom utils packageDescription
#' 
.onAttach <- function(libname, pkgname){
  
  
  desc <- packageDescription("pwrRasch")
  d1 <- desc$Version
  nk <- paste0(rep(" ", 20 - nchar(d1)))
  
  packageStartupMessage("|----------------------------------------------------------|\n",
                          paste0("| ", desc$Package, " ", d1," (",desc$Date,")"), nk, "               |\n" , 
                        "| Statistical Power Simulation for Testing the Rasch Model |\n" ,
                        "|----------------------------------------------------------|" )

}

version <- function(pkg = "pwrRasch") {
  
  lib <- dirname(system.file(package = pkg))
  desc <- packageDescription(pkg)
  
  return(paste(desc$Package, desc$Version, desc$Date,lib))

}
