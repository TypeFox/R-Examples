# THIS FILE CONTAINS PACKAGE HOOKS FOR RNRFA
# ------------------------------------------

# @...: nothing
# spill-over: output information about RNRFA
.onAttach <- function(...) {

  package.name <- "rnrfa"
  mylib <- dirname(system.file(package = package.name))
  ver <- packageDescription(package.name, lib.loc = mylib)$Version
  build.date <- packageDescription(package.name, lib.loc = mylib)$Date


  # build info
  packageStartupMessage("RNRFA (Versions ", ver, ", built: ", build.date, ")")

  # cat, for readability of the message text

  # RNRFA info - do not exceed 80char/line
  packageStartupMessage("
+----------------------------------------------------------------+
|  If you wish to use NRFA Data, please refer to the following   |
|  Terms & Conditions:                                           |
|  http://nrfa.ceh.ac.uk/costs-terms-and-conditions              |
+----------------------------------------------------------------+

")
}
