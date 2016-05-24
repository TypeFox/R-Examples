### Do not delete this file and file name!!
### This file should be loaded before all other *.r files.

### This is to avoid the false positive messages from R CMD check.
###   "no visible binding for global variable"
### Suggested by Prof Brian Ripley
### ?globalVariables

if(getRversion() >= "2.15.1"){
  utils::globalVariables(c(".CF.CT", ".CF.CONF", ".CF.OP", ".CF.DP", ".CF.AC",
                           ".CF.PARAM", ".CO.CT", ".cubfitsEnv",
                           ".CF.GV", ".CF.PT"))
}
