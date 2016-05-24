setMethodS3("rsource", "default", function(file, path=NULL, envir=parent.frame(), output="", buffered=FALSE, ...) {
  rcat(file=file, path=path, envir=envir, output=output, buffered=buffered, ...);
}) # rsource()

setMethodS3("rsource", "RspString", function(..., envir=parent.frame(), output="", buffered=FALSE) {
  rcat(..., envir=envir, output=output, buffered=buffered);
}, protected=TRUE) # rsource()

setMethodS3("rsource", "RspDocument", rsource.RspString)
setMethodS3("rsource", "RspRSourceCode", rsource.RspString)
setMethodS3("rsource", "function", rsource.RspString)
setMethodS3("rsource", "expression", rsource.RspString)

############################################################################
# HISTORY:
# 2014-01-02
# o CLEANUP: Now rstring() methods for several classes uses the exact
#   same function definition.  Also harmonized the ordering of arguments.
# 2013-11-23
# o BUG FIX: rsource() would not evaluate in the current environment.
# 2013-08-05
# o Created.
############################################################################
