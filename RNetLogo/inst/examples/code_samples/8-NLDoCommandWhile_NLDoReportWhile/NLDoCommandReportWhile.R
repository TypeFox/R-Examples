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
model.path <- "/models/Sample Models/Earth Science/Fire.nlogo"
NLLoadModel(paste(path.to.NetLogo,model.path,sep=""), nl.obj=nl.test1)

# submit a single command
NLCommand("setup", nl.obj=nl.test1)

# execute "go" while burned-trees are lower 400 trees and ticks (simulation steps) are lower 30
NLDoCommandWhile("burned-trees < 400 and ticks < 30", "go", nl.obj=nl.test1)

# same as NLDoCommandWhile but save the number of burned-trees and embers (as list with nested lists)
burned.embers <- NLDoReportWhile("burned-trees < 1500 and ticks < 50", "go", c("ticks","burned-trees","count embers"), nl.obj=nl.test1)
print(burned.embers)

# same as before, but save result as data.frame with column names
burned.embers <- NLDoReportWhile("burned-trees < 2500 and ticks < 100", "go", c("ticks","burned-trees","count embers"), as.data.frame=TRUE, df.col.names=c("ticks","burned","embers"), nl.obj=nl.test1)
print(burned.embers)


# use NLQuit(nl.obj=nl.test1) to close the NetLogo Window
