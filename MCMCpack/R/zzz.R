##########################################################################
## start-up and clean-up functions
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################

.onAttach <- function(...) {
 
   # figure out year automatically (probably could be done more elegantly)
   date <- date()
   x <- regexpr("[0-9]{4}", date)
   this.year <- substr(date, x[1], x[1] + attr(x, "match.length") - 1)
   
   # echo output to screen
   packageStartupMessage("##\n## Markov Chain Monte Carlo Package (MCMCpack)")
   packageStartupMessage("## Copyright (C) 2003-", this.year,
      " Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park", sep="")
   packageStartupMessage("##\n## Support provided by the U.S. National Science Foundation")
   packageStartupMessage("## (Grants SES-0350646 and SES-0350613)\n##")
   ## require(coda, quietly=TRUE)
   ## require(MASS, quietly=TRUE)
   ## require(stats, quietly=TRUE)
}

.onUnload <- function(libpath) {
    library.dynam.unload("MCMCpack", libpath)
}

