
##==============================================================================
## Least distance programming
##==============================================================================

Ldei <- function(...) UseMethod ("Ldei")
Ldei.character <- function(...) Ldei.limfile(...)
Ldei.double <- function(...) ldei(...)

##==============================================================================
# reads an input file and solves the least distance model
##==============================================================================

Ldei.limfile <- function(file, verbose=TRUE,...) {
  lim    <- Setup.limfile(file, verbose=verbose)
  Ldei.lim(lim,...)
}

##==============================================================================
## Solves inverse model, least distance programming
##==============================================================================

Ldei.lim <- function(lim,...) {
  ld<-ldei (E=lim$A,F=lim$B,G=lim$G,H=lim$H,...)
  names(ld$X) <- lim$Unknowns
  if (ld$IsError)
     warning("Problem could not be solved")

  return(ld)
}

