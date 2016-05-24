.onLoad <- function(lib, pkg) {
       .C('setNumSlaves', n = as.integer(1L), PACKAGE = "spMC")
}

.onAttach <- function(lib, pkg) {
       if (.Call("isOmp", PACKAGE = "spMC")) {
           packageStartupMessage("Package built with openMP.")
       }
       else {
           packageStartupMessage("Package built without openMP.")
       }
       nn <- setCores()
       if (nn > 0L) packageStartupMessage("Use the function setCores() to change the number of CPU cores.")
}
