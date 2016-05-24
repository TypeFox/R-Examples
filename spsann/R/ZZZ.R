.onAttach <- function (lib, pkg) {
  
  pkg.info <- drop(
    read.dcf(file = system.file("DESCRIPTION", package = "spsann"), fields = c("Title", "Version","Date")))
  
  packageStartupMessage(
    paste("---------------------------------------------------------------\n",
          pkg.info["Title"],                                            " \n",
          "spsann version ", pkg.info["Version"],                       " \n",
          "(built on ", pkg.info["Date"], ") is now loaded                \n",
          "ATTENTION: THIS VERSION CONTAINS BACKWARD INCOMPATIBLE CHANGES!\n",
          "Check the package NEWS and read the new documentation          \n",
          "before you start using the package.                            \n",
          "---------------------------------------------------------------\n",
          sep = "")
  )
}
