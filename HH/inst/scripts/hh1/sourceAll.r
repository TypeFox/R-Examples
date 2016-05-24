## This file described use of the files for the First Edition

## This function is from ?source in R-2.15.0

## If you want to source() a bunch of files, something like
## the following may be useful:
 sourceDir <- function(path, trace = TRUE, ...) {
   ## for (nm in list.files(path, pattern = "\*\\.[RrSsQq]$")) { ## original
    for (nm in list.files(path, pattern = "Ch[01]*\\.[RrSsQq]$")) { ## added "Ch[01].*"
       if(trace) cat(nm,":")
       source(file.path(path, nm), ...)
       if(trace) cat("\n")
    }
 }

## this line sources all r files into R for the 18 chapters of the HH First Edition book
## sourceDir("c:/HOME/rmh/HH-R.package/HH/inst/scripts/hh1", chdir=TRUE, echo=TRUE) ## not working yet


## this line sources all r files into R for the 18 chapters of the HH book
path <- "c:/HOME/rmh/HH-R.package/HH/inst/scripts/hh1"  ## from source directory
path <- system.file("scripts/hh1", package="HH")  ## from installed directory
for (i in list.files(path, pattern = "^Ch[01].*\\.[RrSsQq]$", full.names=TRUE))  ## seems to work
  source(i, chdir=TRUE, echo=TRUE)

## R
## these lines source all r files into R for the 18 chapters of the HH book
##
##  old.HH.ROOT.DIR <- options(HH.ROOT.DIR="c:/HOME/rmh/HH-R.package/HH/inst/scripts/hh1")
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/hh1/Ch02-data.r"                 , chdir=TRUE, echo=TRUE, print.eval=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/hh1/Ch03-conc.r"                 , chdir=TRUE, echo=TRUE, print.eval=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/hh1/Ch04-grap.r"                 , chdir=TRUE, echo=TRUE, print.eval=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/hh1/Ch05-iinf.r"                 , chdir=TRUE, echo=TRUE, print.eval=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/hh1/Ch06-oway.r"                 , chdir=TRUE, echo=TRUE, print.eval=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/hh1/Ch07-mcomp.r"                , chdir=TRUE, echo=TRUE, print.eval=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/hh1/Ch08-rega.r"                 , chdir=TRUE, echo=TRUE, print.eval=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/hh1/Ch09-regb.r"                 , chdir=TRUE, echo=TRUE, print.eval=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/hh1/Ch10-regbb.r"                , chdir=TRUE, echo=TRUE, print.eval=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/hh1/Ch11-regc.r"                 , chdir=TRUE, echo=TRUE, print.eval=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/hh1/Ch12-rhiz-bwplot-alternate.r", chdir=TRUE, echo=TRUE, print.eval=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/hh1/Ch12-tway.r"                 , chdir=TRUE, echo=TRUE, print.eval=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/hh1/Ch13-dsgn.r"                 , chdir=TRUE, echo=TRUE, print.eval=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/hh1/Ch14-apple3.r"               , chdir=TRUE, echo=TRUE, print.eval=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/hh1/Ch14-dsgntwo.r"              , chdir=TRUE, echo=TRUE, print.eval=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/hh1/Ch15-twtb.r"                 , chdir=TRUE, echo=TRUE, print.eval=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/hh1/Ch16-npar.r"                 , chdir=TRUE, echo=TRUE, print.eval=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/hh1/Ch17-logi.r"                 , chdir=TRUE, echo=TRUE, print.eval=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/hh1/Ch18-tser.r"                 , chdir=TRUE, echo=TRUE, print.eval=TRUE)       ##  options(old.HH.ROOT.DIR)





## S-Plus
## these lines source all r files into S-Plus for the 18 chapters of the HH book
## The older version of HH for S-Plus does not have the hh1 subdirectory
##
##  old.HH.ROOT.DIR <- options(HH.ROOT.DIR="c:/HOME/rmh/HH-R.package/HH/inst/scripts")  ## from source directory
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/Ch02-data.r"                 , echo=TRUE, auto.print=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/Ch03-conc.r"                 , echo=TRUE, auto.print=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/Ch04-grap.r"                 , echo=TRUE, auto.print=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/Ch05-iinf.r"                 , echo=TRUE, auto.print=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/Ch06-oway.r"                 , echo=TRUE, auto.print=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/Ch07-mcomp.r"                , echo=TRUE, auto.print=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/Ch08-rega.r"                 , echo=TRUE, auto.print=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/Ch09-regb.r"                 , echo=TRUE, auto.print=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/Ch10-regbb.r"                , echo=TRUE, auto.print=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/Ch11-regc.r"                 , echo=TRUE, auto.print=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/Ch12-rhiz-bwplot-alternate.r", echo=TRUE, auto.print=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/Ch12-tway.r"                 , echo=TRUE, auto.print=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/Ch13-dsgn.r"                 , echo=TRUE, auto.print=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/Ch14-apple3.r"               , echo=TRUE, auto.print=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/Ch14-dsgntwo.r"              , echo=TRUE, auto.print=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/Ch15-twtb.r"                 , echo=TRUE, auto.print=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/Ch16-npar.r"                 , echo=TRUE, auto.print=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/Ch17-logi.r"                 , echo=TRUE, auto.print=TRUE)
source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/Ch18-tser.r"                 , echo=TRUE, auto.print=TRUE)
##  options(old.HH.ROOT.DIR)
##  rm(old.HH.ROOT.DIR)
