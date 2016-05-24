################################################################################
#
# TAQ - 
#
# collected by Fabrizio Cipollini
# Version 
#
################################################################################

.onLoad <-
function(lib, pkg)
{
  ## FUNCTION:

  #### Load dll:
  library.dynam(chname = .package.name(), package = pkg, lib.loc = lib)
  
  #### Answer
  invisible(NULL)
}
# ------------------------------------------------------------------------------


.onAttach <-
function(lib, pkg)
{
  ## FUNCTION:

  #### Package:
  packageStartupMessage( 
    "---------------------------------------------------------------------\n", 
    paste(
      "--                             ", .package.name(), 
      "                             --\n", sep = ""),  
      "---------------------------------------------------------------------\n", 
    "    Package attached\n",
    domain = NULL, appendLF = TRUE)

  #### Load machine constants
  # .Machine.constants()

  #### Answer
  invisible(NULL)
}
# ------------------------------------------------------------------------------


.onUnload <-
function(libpath = NULL)
{
  ## FUNCTION:
  
  #### Unload dll:
  library.dynam.unload(chname = .package.name(), libpath)

  #### Answer
  invisible(NULL)
}
# ------------------------------------------------------------------------------


.Last.lib <-
function(libpath = NULL)
{
  ## FUNCTION:
  
  #### Package:
  cat("---------------------------------------------------------------------\n")
  cat("--                           ", .package.name(), 
    "                                   --\n")
  cat("---------------------------------------------------------------------\n")
  cat("                      Package detached\n")
}
# ------------------------------------------------------------------------------
