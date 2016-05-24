checklogicalbounds <- function(v) { ## just check logical bounds on the kriged variables.
  varnames <- colnames(blackbox.getOption("fitobject")$x)
  for (vit in varnames) {
    blob <- v[vit]
    if(islogscale(vit)) blob <- exp(blob)
    if (vit=="twoNmu") {if (blob<=0) {return(FALSE);} else {next;}}
    if (vit=="twoNm") {if (blob<=0) {return(FALSE);} else {next;}}
    if (vit=="g") {if (blob<0 || blob>1) {return(FALSE);} else {next;}}
    if (vit=="latt2Ns2") {if (blob<=0) {return(FALSE);} else {next;}}
    ## if we reach this point then...
    stop.redef(paste("(!) Unhandled variable", vit, "in checklogicalbounds()"))
  }
  ## if we reach this point then...
  return(TRUE)
} ## end def checklogicalbounds
