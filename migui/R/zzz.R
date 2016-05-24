.onLoad <- function(lib, pkg) {
  return(invisible(NULL))
}

.onUnLoad <- function(lib, pkg) {
  return(invisible(NULL))
}

.onAttach <- function(...) {
  miguiLib <- dirname(system.file(package = "migui"))
  version <- utils::packageDescription("migui", lib.loc = miguiLib)$Version
  packageStartupMessage(paste("\nmigui (Version ", version, ")", sep = ""))
  packageStartupMessage("migui Copyright (C) 2010, 2011, Columbia University")
  packageStartupMessage("This program comes with ABSOLUTELY NO WARRANTY.")
  packageStartupMessage("This is free software, and you are welcome to redistribute it")
  packageStartupMessage("under the General Public License version 2 or later.")
  packageStartupMessage("Execute RShowDoc('COPYING') for details.")
  return(invisible(NULL))
}
