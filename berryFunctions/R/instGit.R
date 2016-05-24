#' install github package
#'
#' Quickly install a package from github without having to install devtools with all its dependencies.
#'
#' @details Works only for pure R package structure repositories from the master branch.
#' Installs package dependencies listed in 'Imports' and 'Depends', but ignores version requirements!
#' Tested only on windows 7 with R3.2.2. 
#' Note: \code{devtools::install_github} is much more extensive! \cr
#' Note: drat is also much better than this quick hack.
#' \url{http://dirk.eddelbuettel.com/code/drat.html},
#' \url{https://github.com/eddelbuettel/drat},
#' \url{http://eddelbuettel.github.io/drat/DratForPackageAuthors.html}
#' Give your github users this code:\cr
#' \code{source("https://raw.githubusercontent.com/brry/berryFunctions/master/R/instGit.R")}\cr
#' \code{instGit("brry/extremeStat")}\cr
#' \code{library(extremeStat)}\cr
#'
#' @return NULL
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Dec 2015 + Mar/Apr 2016
#' @seealso \code{\link{funSource}}, \code{install_github} in each of the packages \code{devtools, ghit, remotes}
#' @keywords package
#' @export
#' @examples
#'
# Worked fine on my computer for:
#' if(FALSE){
#' instGit("talgalili/installr")
#' instGit("talgalili/installr", FALSE)
#' instGit("hadley/readxl")
#' instGit("mages/googleVis") # many dependencies!
#' instGit("twitter/AnomalyDetection")
#' instGit("yihui/knitr")
#' instGit("ramnathv/slidify")
#' instGit("jrnold/ggthemes")
#' }
#'
#' @param pk Character string in the from of "user/package"
#' @param cleanup Remove downloaded zipfile and folder with source code. DEFAULT: TRUE
#' @param \dots Further arguments passed to \code{\link{install.packages}}, untested so far
#'
instGit <- function(
  pk,
  cleanup=TRUE,
  ...)
{
owd <- setwd(tempdir())
on.exit(setwd(owd), add=TRUE)
pkn <- strsplit(pk, "/")[[1]][2] # package name part
# Download the zip file with the source code:
htt <- if(getRversion() < "3.2.2") "http" else "https"
suppressWarnings(
download.file(url=paste0(htt,"://github.com/",pk,"/archive/master.zip"),
              destfile=paste0(pkn,".zip"))
) # suppress warnings like downloaded length 350226 != reported length 200
# unzip and rename the folder:
unzip(paste0(pkn,".zip"))
file.rename(paste0(pkn,"-master"), pkn)
# Find out dependencies - really not elegant at all!:
deps <- read.dcf(paste0(pkn, "/DESCRIPTION"), fields="Imports")
deps2<- read.dcf(paste0(pkn, "/DESCRIPTION"), fields="Depends")
deps <- paste(deps,deps2, sep=",")
deps <- gsub("\n", "", deps) # remove line breaks
while(grepl(" ", deps)) deps <- gsub(" ", "", deps) # remove spaces
deps <- gsub("NA", "", deps) # remove NAs for packages not listing fields
deps <- strsplit(deps, ",")[[1]] # split entries
deps <- deps[deps!=""] # NA leftover
deps <- sapply(strsplit(deps, "(", fixed=T), "[", 1) # remove version restrictions
deps <- deps[deps!="R"] # remove R (often in depends with version)
isinst <- deps %in% rownames(installed.packages())
depsinst <- deps[isinst]
deps <- deps[!isinst] # install only new packages
# tell user about installing dependencies:
if(!all(is.na(depsinst)))
   message("--- The following dependencies were already installed, but versions are unchecked: \n ",
           paste(depsinst, collapse=", "))
if(!all(is.na(deps)))
   message("--- instGit will now install the following dependencies: \n ",
           paste(deps, collapse=", "))
flush.console()
# install dependencies:
dummy <- lapply(deps, install.packages, ...)
# actually install the package itself:
message("--- instGit will now install ", pkn, " ...")
flush.console()
install.packages(pkn, repos=NULL, type="source", ...)
if(!cleanup) message("--- The downloaded zip is in: ", getwd())
# clean up:
if(cleanup)
  {
  unlink(pkn, recursive=TRUE)
  unlink(paste0(pkn,".zip"))
  }
} #end of function
