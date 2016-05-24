MRSPwelcomemessage <- function(){
  cat("-----------------------------",
      "This is MRSP version 0.4.3",
      #"-----------------------------",
      "Author: Wolfgang Poessnecker",
      "Date: December 09, 2014",
      "-----------------------------",
      sep = "\n")   }

.onAttach <- function(libname, pkgname){
  packageStartupMessage(MRSPwelcomemessage())
}
