#
# This script runs all unit tests in this package
#
# set the working directory to the location of R source code 
if (!is.null(getOption('rworkspace', NULL))) {
  systemInfo <- Sys.info() # get system info
  workingdir <- paste(getOption('rworkspace'),"causata/R",sep="/")
  setwd(paste(getOption('rworkspace'),"causata/R",sep="/"))
}

#library(causata) # do not load the causata library as it may be out of sync with the latest code in the files
library(testthat)

# clear memory before testing
rm(list=ls())

# source the files to be tested so that we are testing the latest versions, not necessarily in sync
# with the library
#for (file in dir()) source(file) 
  
# run tests
test_dir('tests', reporter='Summary')