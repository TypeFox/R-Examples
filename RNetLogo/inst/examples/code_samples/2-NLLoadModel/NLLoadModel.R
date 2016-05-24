library(RNetLogo)

# PLEASE FILL IN THE PATH TO YOUR NetLogo INSTALLATION FOLDER !!!!
path.to.NetLogo <- "C:/Program Files/NetLogo 5.3/app"

# start NetLogo, if it is not started yet
if (!exists("nl.test1", -1)) 
{
  # an reference name
  nl.test1 <- "nl.test1"
  NLStart(path.to.NetLogo, gui=TRUE, nl.obj=nl.test1)
}

model.path <- "/models/Sample Models/Earth Science/Fire.nlogo"
NLLoadModel(paste(path.to.NetLogo,model.path,sep=""),nl.obj=nl.test1)


# use NLQuit(nl.obj=nl.test1) to close the NetLogo Window
