## Jonathan Wand : wand(at)stanford.edu
## a little helper function for demos
allequal.test <- function( target, current, expect=TRUE) {
  z <- all.equal(target,current)
  if (typeof(z) != "logical")
    z <- FALSE

  if (expect == z) { 
    cat("Expecting all.equal() =",expect," :  test is correctly",z,"\n")
  } else {
    stop(paste("Was EXPECTING test to be",expect," but all.equal() =",z,"\n"))
  }
  return(invisible(z))
}
