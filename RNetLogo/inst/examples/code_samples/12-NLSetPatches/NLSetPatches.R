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

# get the world dimension
world.dim <- NLReport(c("world-width","world-height"), nl.obj=nl.test1)
print(world.dim)

# create a vector with length: world-width * world-height and incremented values
my.vector <- 1:(world.dim[[1]]*world.dim[[2]])
print(my.vector)

# create a matrix with number of columns equiv. to world-width
my.matrix <- matrix(my.vector, ncol=world.dim[[1]])
print(my.matrix)

# use this matrix as input for "plabel"
NLSetPatches("plabel", my.matrix, nl.obj=nl.test1)

             
# use NLQuit(nl.obj=nl.test1) to close the NetLgogo Window
