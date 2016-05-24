.onAttach <- function(...) {
  welcome <- paste(""                                              ,
                   "----------------------------------------------",
                   "  'blockseg' package version 0.1-0           ",
                   ""                                              ,
                   " Still under development... feedback welcome  ",
                   "----------------------------------------------",
                   sep = "\n")
  packageStartupMessage(welcome)
}

