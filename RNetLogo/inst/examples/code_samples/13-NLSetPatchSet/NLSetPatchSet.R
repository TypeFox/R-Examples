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
path.to.local.file <- 'C:/Users/jthiele/Documents/R/win-library/3.2/RNetLogo/examples/code_samples/12-NLSetPatches/patchtest.nlogo'
NLLoadModel(path.to.local.file, nl.obj=nl.test1)
                    
# submit a single command
NLCommand("setup", nl.obj=nl.test1)

# get some patches with their xpcor, ypcor and pcolor
some.patches <- NLGetPatches(c("pxcor","pycor","pcolor"), "patches with [plabel > 15]", nl.obj=nl.test1)

# change the value of pcolor to a random value
some.patches$pcolor <- floor(abs(rnorm(nrow(some.patches))*100)) 

# send the updated pcolor values for the selected pathces 
# to NetLogo 
patch.var <- "pcolor"
NLSetPatchSet(some.patches, patch.var, nl.obj=nl.test1)
             
# use NLQuit(nl.obj=nl.test1) to close the NetLgogo Window
