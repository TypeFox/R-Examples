provideDevice <- function(file=blackbox.getOption("basicRplotsfile"),bbdefaultPars,...) {
  if (is.null(file)) {
    ## typically the case in interactive session unless the user (or previous code) defined a basicRplotsfile in such a session
    dev <- getOption("device")
    if (  ! (class(dev)=="character" && dev == "RStudioGD") ) dev.new(...)
    if (bbdefaultPars) par(blackbox.getOption("ParArgs"))
  } else {
    newFile <- providePlotFile(file) ## will check if file is not already opened and, if not,
    ## providePlotFile    will eval(call(blackbox.getOption("graphicsFormat"), file,...))
    if (newFile && bbdefaultPars) par(blackbox.getOption("ParArgs"))
  }
}
