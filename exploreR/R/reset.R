#' Reset R
#'
#' This function simply erases the the console, detaches all packages and removes all data from the global environment.  The purpose of which is to provide an easy command which can clean up the workspace.  Very useful after you spend some time experimenting.
#'
#' @export
#' @examples 
#' reset()

reset <- function() {
  
  # This function is exactly the same as rm(list = ls()), but because it happens within a function, the environment has to be declared as the global environment.
  rm(list = ls(envir = globalenv()), envir = globalenv())
  
  # This call sends the Ctrl + L function to the console.
  cat("\014")
  
  # Generate the list of currently installed packages.
  installed.packages <- .packages()
  
  # This vector is all the packagges I would like to leave installed.  "simpletools" in the name of the package that contains this function, the rest are default packages and included to remove the possbility of generating unforeseen errors.
  safe.packages <- c("stats", "graphics", "grDevices", "utils", "datasets", "methods", "base", "exploreR")
  
  # For loop that goes through every package installed.
  for (i in installed.packages) {
    
    # If the package DOES NOT appear on the safe list, kill it.
    if (!i %in% safe.packages) { 
      
      kill.package <- paste("package", i, sep = ":")
        
        detach(kill.package, unload = TRUE, character.only = TRUE)
   }
  
  }
  
}