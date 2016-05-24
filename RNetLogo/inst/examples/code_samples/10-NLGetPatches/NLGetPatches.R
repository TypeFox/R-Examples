library(RNetLogo)

# path to NetLogo installation folder
# PLEASE FILL IN THE PATH TO YOUR NetLogo INSTALLATION FOLDER !!!!
path.to.NetLogo <- "C:/Program Files/NetLogo 5.3/app"

# start NetLogo, if it is not started yet
if (!exists("nl.test1", -1)) 
{
  # an reference name
  nl.test1 <- "nl.test1"
  NLStart(path.to.NetLogo, gui=TRUE, nl.obj=nl.test1)
}

# load a sample model
# PLEASE FILL IN THE PATH TO THE SAMPLE NetLogo FILE !!!!
path.to.local.file <- 'C:/Users/jthiele/Documents/R/win-library/3.2/RNetLogo/examples/code_samples/10-NLGetPatches/patchtest.nlogo'
NLLoadModel(path.to.local.file, nl.obj=nl.test1)

# submit a single command
NLCommand("setup", nl.obj=nl.test1)

# get the values of all patches (the world) of the variable plabel as data.frame
# orderd from upper left to bottom right
plabels <- NLGetPatches("plabel", nl.obj=nl.test1)
print(plabels)

# get the same values, but as matrix, looking like the NetLogo world 
# (but maybe with different indices, because R matrix always starting with (1,1) in upper left)
plabels.matrix <- NLGetPatches("plabel", as.matrix=TRUE, nl.obj=nl.test1)
print(plabels.matrix)

# get a subset of patches, e.g. all patches with values of plabel greater than 10 together with the x and y coordinates 
plabels.gt.10 <- NLGetPatches(c("plabel","pxcor", "pycor"), "patches with [plabel > 10]", nl.obj=nl.test1)
print(plabels.gt.10)

# same as before but as list (where each list element represent a patch variable)
plabels.gt.10.list <- NLGetPatches(c("plabel","pxcor", "pycor"), "patches with [plabel > 10]", as.data.frame=FALSE, nl.obj=nl.test1)
print(plabels.gt.10.list)

# same as before but as a list where each list element represent a patch 
# (ATTENTION: this is very slow! Use it carefully, especially on a large number of patches.)
plabels.gt.10.list2 <- NLGetPatches(c("plabel","pxcor", "pycor"), "patches with [plabel > 10]", as.data.frame=FALSE, patches.by.row=TRUE, nl.obj=nl.test1)
print(plabels.gt.10.list2)

# use NLQuit(nl.obj=nl.test1) to close the NetLogo Window
