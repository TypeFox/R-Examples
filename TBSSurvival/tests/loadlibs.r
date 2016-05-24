# TBSSurvival package for R (http://www.R-project.org)
# Copyright (C) 2013 Adriano Polpo, Cassio de Campos, Debajyoti Sinha
#                    Jianchang Lin and Stuart Lipsitz.
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

## This code is used for testing purposes. The TBSSurvival library does not
## depend on it for any of its functionalities

installpacks <- function(loc=NULL,repos="http://stat.ethz.ch/CRAN/") {
  ## set the repository to use
  options(repos=repos)
  ## install the packages
  install.packages("coda",lib=loc)
  install.packages("mcmc",lib=loc)
  install.packages("normalp",lib=loc)
  install.packages("R.methodsS3",lib=loc)
  install.packages("R.oo",lib=loc)
  install.packages("R.utils",lib=loc)
  install.packages("Rsolnp",lib=loc)
  install.packages("survival",lib=loc)
#  install.packages("e1071",lib=loc)
#  install.packages("eha",lib=loc)
  install.packages("truncnorm",lib=loc)
  install.packages("BMS",lib=loc)
  
  ## this following line install the TBS package itself, so nothing else is needed.
  ## For testing, sometimes it is better to work without installing it for a while...
  ##      install.packages('./TBSSurvival_version.tar.gz',repos=NULL,type="source")
}

loadlibs <- function(libdir=NULL) {
  w <- options("warn")
  options("warn" = -1)
  if (require("TBSSurvival",quietly=TRUE)==FALSE) {
    library("BMS",lib.loc=libdir)
    library("coda",lib.loc=libdir)
    library("mcmc",lib.loc=libdir)
    library("normalp",lib.loc=libdir)
    library("R.methodsS3",lib.loc=libdir)
    library("R.oo",lib.loc=libdir)
    library("R.utils",lib.loc=libdir)
    library("Rsolnp",lib.loc=libdir)
    library("survival",lib.loc=libdir)
#    library("e1071",lib.loc=libdir)
#    library("eha",lib.loc=libdir)
    library("truncnorm",lib.loc=libdir)
    source("../R/tbs.survreg.be.r")
    source("../R/ptbs.r")
    source("../R/qtbs.r")
    source("../R/dtbs.r")
    source("../R/rtbs.r")
    source("../R/htbs.r")
    source("../R/tbs.survreg.mle.r")
    source("../R/local.r")
    source("../R/dt2.r")
    source("../R/dlogis2.r")
    source("../R/dist.error.r")
    
  } else {
    library("TBSSurvival")
  }
  options("warn" = w[[1]])
}

## Load data
alloyT7987 <- read.table("../data/alloyT7987.txt",header=TRUE)

