# files and functions that should be loaded at last when calling
# library(LambertW)

.onAttach <- function(libname, pkgname){
  welcome.msg <- 
    paste0("This is 'LambertW' version ", utils::packageVersion("LambertW"), 
           '.  Please see the NEWS file and citation("LambertW").\n')
  packageStartupMessage(welcome.msg, domain = NULL, 
                        appendLF = TRUE)
}