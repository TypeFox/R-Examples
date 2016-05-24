### R code from vignette source 'range_expansion.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: range_expansion.Rnw:29-43
###################################################
library(MetaLandSim)

#Load starting landscape (the simulation will assume that 
#all subsequent landscapes are built with the same parameter combination).

data(rland)

#Create range expansion model. Here run only with two repetitions (iter=2). 
#Ideally it should be run with more repetitions to provide more robust results.

data(param1)

rg_exp1 <- range_expansion(rl=rland, percI=50,  param=param1,b=1, tsteps=100, 
iter=2)


###################################################
### code chunk number 2: range_expansion.Rnw:59-78 (eval = FALSE)
###################################################
## data(rg_exp)
## presences <- paste(system.file(package="MetaLandSim"),
##  "/examples/presences.asc", sep="")
## landmask <- paste(system.file(package="MetaLandSim"), 
## "/examples/landmask.asc", sep="")
## 
## library(rgrass7)
## 
## #First, start GRASS from R: 
## initGRASS(gisBase = "grass folder", home = tempdir(), 
## gisDbase = "mapset location",override = TRUE)
## 
## #Create raster, using the sample dataset 
## #rg_exp (generated with 100 repetitions)
## 
## data("rg_exp")
## 
## range_raster(presences.map = presences, re.out=rg_exp, 
## mask.map=landmask, plot.directions=FALSE)


