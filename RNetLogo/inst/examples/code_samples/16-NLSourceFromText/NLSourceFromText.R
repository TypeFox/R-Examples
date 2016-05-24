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

# with NLSourceFromText we can create or append NetLogo Model source from R.
# This function is only available when NetLogo was started with GUI (not in headless mode)

# load a sample model
# PLEASE FILL IN THE PATH TO THE SAMPLE NetLogo FILE !!!!
path.to.local.file <- 'C:/Users/jthiele/Documents/R/win-library/3.2/RNetLogo/examples/code_samples/16-NLSourceFromText/samplemodel.nlogo'
NLLoadModel(path.to.local.file, nl.obj=nl.test1)

# append the model by a new procedure called "go2" (there is no way to change the current model)
# which creates a new turtle and calls go
go2 <- "to go2\n crt 1\n go\nend"
NLSourceFromString(go2, nl.obj=nl.test1)

# initialize the model
NLCommand("setup", nl.obj=nl.test1)

# run the model and save the number of turtles (to check, if our new procedure works correct)
no.turtles <- NLDoReport(10, "go2", c("ticks","count turtles"), as.data.frame=TRUE, df.col.names=c("time","turtles"), nl.obj=nl.test1)
plot(no.turtles, type='l')


# the second way is, to create a completely new model
setup <- "to setup\n ca\n reset-ticks\nend"
go <- "to go\n ask patches [\n  set pcolor random 255\n ]\n tick\nend"

NLSourceFromString(c(setup,go), append.model=FALSE, nl.obj=nl.test1)

# initialize and run the new model
NLCommand("setup", nl.obj=nl.test1)
NLDoCommand(1000, "go", nl.obj=nl.test1)

# we can now change the go procedure, submit the source again and can continue our simulation without resetting it
# with this method you can create NetLogo model code depending on some other calculations in R
go <- "to go\n ask patches [\n  set plabel random 255\n ]\n tick\nend"
NLSourceFromString(c(setup,go), append.model=FALSE, nl.obj=nl.test1)
NLDoCommand(1000, "go", nl.obj=nl.test1)

# use NLQuit(nl.obj=nl.test1) to close the NetLogo Window
