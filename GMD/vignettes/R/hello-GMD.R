##' Check GMD's sanity and start up
##' @name hello-GMD

## GMD at CRAN, for source code download and installation
## http://cran.r-project.org/web/packages/GMD/index.html

## load GMD
library(GMD)

## version of GMD and description
packageVersion("GMD")
packageDescription("GMD")

## view GMD vignette
vignette("GMD-vignette",package="GMD")

## list the available data sets in GMD
data(package="GMD")

## list all the objects in the GMD
ls("package:GMD")

## help info on GMD
help(package="GMD")

## run a demo
demo("GMD-demo")

## cite GMD in publications
citation(package="GMD")

