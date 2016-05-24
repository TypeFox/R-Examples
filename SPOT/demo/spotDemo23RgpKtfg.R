##
## Use demo(... , ask=F) to run the demo all at once, without interruption.
##
## This example shows how to tune RGP with a Multiple Algorithm Multiple
## Problems (MAMP) approach and mixed model analysis. A Kriging-based test
## function generator (KTFG) is used to create randomized variants of a
## symbolic regression test function. These instances are then solved 
## with different combinations of RGP parameter settings. Finally, a mixed
## model analysis is performed on the results.
##

## Load SPOT package...
require("SPOT")

## Load rgp
require("rgp")

## Get path of test project...
spotPath <- find.package("SPOT")
demoPath <- file.path(spotPath, "demo23RgpKtfg")

source(file.path(demoPath,"spotAlgStartKtfgRgp.R"))

## Show path...
message(demoPath)

## Get the SPOT config file...
confFile <- file.path(demoPath, "ktfgRgp001.conf")

## Initialize the design to be evaluated...
res <- spot(confFile, spotTask = "init")

## Run algorithm (RGP) on KTFG-generated problem instance...
## please note that this might take several seconds
res <- spot(confFile, spotTask = "run", spotConfig = res)

## Generate a report...
res <- spot(confFile, spotTask = "rep", spotConfig = res)

## EOF 

