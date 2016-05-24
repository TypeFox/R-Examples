### Do not delete this file and file name!!
### This file should be loaded before all other *.r files.

### This is to avoid the fale positive messages from R CMD check.
###   "no visible binding for global variable"
### Suggested by Prof Brian Ripley
### ?globalVariables

if(getRversion() >= "2.15.1"){
  utils::globalVariables(c(
    ".EMC",
    ".EMC.Rnd",
    "col.ppcontour",
    "color.class"
    ))
}
