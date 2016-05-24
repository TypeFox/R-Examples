.onAttach <- function(lib,pkg)
{


   # figure out year automatically (probably could be done more elegantly)
   date <- date()
   x <- regexpr("[0-9]{4}", date)
   this.year <- substr(date, x[1], x[1] + attr(x, "match.length") - 1)

   # echo output to screen
   packageStartupMessage("##\n## On Bayesian analysis of Threshold autoregressive model (BAYSTAR)\n")
   packageStartupMessage("## Copyright (C) 2007-", this.year,
      " Cathy W. S. Chen, Edward M.H. Lin, F.C. Liu, and Richard Gerlach\n", sep="")
ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
     ver <- as.character(ver)	
    packageStartupMessage("Version: ", ver, "loaded\n")

#   require(mvtnorm, quietly=TRUE)
#   require(coda, quietly=TRUE)
}