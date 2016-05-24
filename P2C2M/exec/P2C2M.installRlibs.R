#!/usr/bin/R
#author = "Michael Gruenstaeudl, PhD"
#copyright = "Copyright (C) 2014 Michael Gruenstaeudl"
#contributors = c("Michael Gruenstaeudl","Noah Reid")
#email = "gruenstaeudl.1@osu.edu"
#version = "2014.10.06.0200"

p2c2m.install = function() {
  # Descr:  install thir-party executables and libraries necessary for p2c2m
  # Deps:   p2c2m.pkgload
  # I/p:    (none)

 ## Prepare libraries to be loaded
  dwnlList = list()
  # Update the download list periodically (Last update: 2014.09.17)
  dwnlList[["genealogicalSorting"]] = 
    "http://www.genealogicalsorting.org/resources/files/genealogicalSorting_0.92.tar.gz"
  dwnlList[["phybase"]] = 
    "http://odyssey.bioinformatics.uga.edu/~lliu/phybase/phybase_1.3.1.tar.gz"

  dependencies = c("ape", "apTreeshape", "genealogicalSorting", "ggplot2", 
                   "Rmpi", "rPython", "phybase", "stringr", "xtermStyle",
                   "P2C2M")

 ## Set a default CRAN mirror
  options(repos = c(CRAN = "http://cran.rstudio.com"))
  # Other potential selections:
  # http://cran.case.edu/ (in Ohio)
  # http://cran.at.r-project.org (main CRAN mirror)

 ## Load the libraries
  lapply(dependencies, p2c2m.pkgload, dwnlList)

 ## Complete notice
    cat("\n -- P2C2M: Installation of R packages complete --\n\n")
}


p2c2m.pkgload = function(pkgName, dwnlList) {
  # Descr:  loads and, if not found in local library, install packages
  # Deps:   (none)
  # I/p:    pkgName = name of package to be checked/installed
  #         dwnlList = download list

#################################
# 1. Presence of pkg in library #
#################################
  # Note: "require" has to be executed to see if it cannot be executed
  #       Hence, there is no need to write a second line "require(...)"
  if (!(suppressWarnings(require(pkgName, character.only=T)))) {

##################################
# 2. Confirm internet connection #
##################################
    # If a package is not found in the local repository, it must be installed
    # Hence, an internet connection must exit.
    # Checking if connected to internet.
    tryCatch(
      {msg = system("ping -c3 www.google.com", intern=T)}, 
      warning=function(w) {msg=NULL}
    )
    # If not connected to the internet, abort
    if (is.null(msg)) {
      stop(cat("\n\n -- Please connect to the internet. Aborting. --\n\n"))
    }

######################
# 3. Install package #
######################
    cat(paste("\n -- Package '", pkgName, "' needs to be installed. --\n\n"))
    allow = readline(prompt=" May I install the above package? ([y]/n): ")
    if ("N" %in% toupper(allow) | "NO" %in% toupper(allow)) {
      stop(cat(paste("\n\n -- Please install package", pkgName,
                     "manually. Aborting. --\n\n")))
    }
    else {
      # Check if a package is listed on CRAN
      if(!pkgName %in% available.packages()) {
        download.file(dwnlList[[pkgName]], pkgName)
        install.packages(pkgName, dependencies=T, type="source", repos=NULL)
      }
      else {
        install.packages(pkgName, dependencies=T)
      }
      cat("\n")
    }
  }

###################
# 4. Load package #
###################
  library(pkgName, character.only=T)

}
