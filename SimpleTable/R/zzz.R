##########################################################################
## start-up and clean-up functions
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## Copyright (C) 2008 Kevin M. Quinn
##########################################################################

.onAttach <- function(...) {
 
   # figure out year automatically (probably could be done more elegantly)
   date <- date()
   x <- regexpr("[0-9]{4}", date)
   this.year <- substr(date, x[1], x[1] + attr(x, "match.length") - 1)
   
   # echo output to screen
   packageStartupMessage("\n\n    //\n   // Bayesian Inference and Sensitivity Analysis for Causal Effects \n",
       "  //        from 2 x 2 and 2 x 2 x K Tables in the Presence\n",
       " //             of Unmeasured Confounding (SimpleTable)\n", sep="")
   packageStartupMessage("//                 Copyright (C) 2008", #this.year,
      " Kevin M. Quinn\n\n", sep="")
#   cat("##\n## Support provided by the U.S. National Science Foundation\n")
#   cat("## (Grants SES-0350646 and SES-0350613)\n##\n")
   #require(MCMCpack, quietly=TRUE)
   #require(hdrcde, quietly=TRUE)
   #require(locfit, quietly=TRUE)
   #require(tcltk, quietly=TRUE)
}



