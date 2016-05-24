".onAttach" <-
function(lib, pkg) {
  MSG <- packageStartupMessage # renaming the function
  dsc <- packageDescription(pkg) # function will not be found without
  # loading the utils package in the Depends: section of the DESCRIPTION
  if(interactive() || getOption("verbose")) {
    MSG(sprintf("# Hydraulics package %s (%s) loaded.", pkg, dsc$Version))
  }
}
