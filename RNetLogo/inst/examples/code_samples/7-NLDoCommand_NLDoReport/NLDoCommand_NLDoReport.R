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
NLLoadModel(paste(path.to.NetLogo,model.path,sep=""),nl.obj=nl.test1)

# submit a single command
NLCommand("setup", nl.obj=nl.test1)

# execute "go" command 10 times
NLDoCommand(10, "go", nl.obj=nl.test1)

# reset simulation
NLCommand("setup", nl.obj=nl.test1)

# execute "go" command 10 times and report "burned-trees" after each execution
burned.trees <- NLDoReport(10, "go", "burned-trees", nl.obj=nl.test1)
print(burned.trees)

# get more than one value in each step (results in a list with nested list)
burned.percentage.trees <- NLDoReport(10, "go", c("burned-trees", "burned-trees / initial-trees * 100"), nl.obj=nl.test1)
print(burned.percentage.trees)

# get the result as data.frame
burned.percentage.trees.df <- NLDoReport(10, "go", c("burned-trees", "burned-trees / initial-trees * 100"), as.data.frame=TRUE, nl.obj=nl.test1)
print(burned.percentage.trees.df)

# set columnnames for data.frame during execution
burned.percentage.trees.df <- NLDoReport(10, "go", c("burned-trees", "burned-trees / initial-trees * 100"), as.data.frame=TRUE, df.col.names=c('burned','percent'), nl.obj=nl.test1)
print(burned.percentage.trees.df)



# use NLQuit(nl.obj=nl.test1) to close the NetLogo Window
