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

# load a sample model
model.path <- "/models/Sample Models/Earth Science/Fire.nlogo"
NLLoadModel(paste(path.to.NetLogo,model.path,sep=""),nl.obj=nl.test1)

# setup the model and run one step 
NLCommand("setup", "go", nl.obj=nl.test1)

# get a single value
burned <- NLReport("burned-trees", nl.obj=nl.test1)
print(burned)

# get two/multiple values as list
burned.init <- NLReport(c("burned-trees","initial-trees"), nl.obj=nl.test1)
print(burned.init)
  # or
burned.init <- NLReport("(list burned-trees initial-trees)", nl.obj=nl.test1)
print(burned.init)

# some other examples
no.turtles <- NLReport("count turtles", nl.obj=nl.test1)
print(no.turtles)
no.patches <- NLReport("count patches with [pxcor < 10 and pycor > 4]", nl.obj=nl.test1)
print(no.patches)
world.dims <- NLReport(c("world-width","world-height","min-pxcor","max-pxcor","min-pycor","max-pycor"), nl.obj=nl.test1)
print(world.dims)



# use NLQuit(nl.obj=nl.test1) to close the NetLogo Window
