
# ----------------------------------------------------------------
# $Author: thm $
# $Date: 2014-11-21 10:49:55 +0100 (Fri, 21 Nov 2014) $
# $Rev: 321 $
# ----------------------------------------------------------------

.onAttach <- function(libname, pkgname) {
  ## Action taking place when calling 'library(wux)'
  ## Currently it prompts the WUX version number
  ## This is an internal function.
  ##
  ## History:
  ##   2011-09-21 | orig code (thm)
  ##   2014-11-21 | platform independent (thm)
  ## 
  path <- file.path(libname, pkgname)
  description.file <- list.files(path, pattern = "DESCRIPTION")

  if (length(description.file) != 0) {
    ## if description file...
    description.file.full <- file.path(path, description.file)
    description.file.cont <- readLines(description.file.full)
    version <- grep("Version: ", description.file.cont, value = TRUE)
    version.number <- sub("Version: ", "", version)
    
    packageStartupMessage("\nWUX runtime version: ", version.number)
  }
  ## else no prompt
}
