
# hdi method for class function, where 'object' must be an
#   inverse cumulative density function (ICDF).
# Based on Kruschke (2011) Doing Bayesian Data Analysis.

# The ICDF is checked and an informative error message given if it fails.

# The TeachingDemos package by Greg Snow contains a similar function, but
#   does not deal correctly with credMass < 0.5: it returns a CRI with
#   credMass = 1 - user-specified credMass, with no indication of
#   any issue.

hdi.function <- function(object, credMass=0.95, tol, ...)  {
  checkCredMass(credMass)
  if(missing(tol))
    tol <- 1e-8
  if(class(try(object(0.5, ...), TRUE)) == "try-error")
    stop(paste("Incorrect arguments for the inverse cumulative density function",
        substitute(object)))
  # cf. code in Kruschke 2011 p630
  intervalWidth <- function( lowTailPr , ICDF , credMass , ... ) {
    ICDF( credMass + lowTailPr , ... ) - ICDF( lowTailPr , ... )
  }
  optInfo <- optimize( intervalWidth , c( 0 , 1.0 - credMass) , ICDF=object ,
                      credMass=credMass , tol=tol , ... )
  HDIlowTailPr <- optInfo$minimum
  result <- c(lower = object( HDIlowTailPr , ... ) ,
            upper = object( credMass + HDIlowTailPr , ... ) )
  attr(result, "credMass") <- credMass
  return(result)
}
