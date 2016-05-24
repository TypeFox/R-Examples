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
NLLoadModel(paste(path.to.NetLogo,model.path,sep=""), nl.obj=nl.test1)

# submit a single command
NLCommand("setup", nl.obj=nl.test1)

# submit two/multiple commands at once
NLCommand("set density 60", "go", nl.obj=nl.test1)

# submit a vector (overrides the variable burned-trees - just for testing, makes no sense)
my.vector <- c(1,2,3,4,5,6)
NLCommand("set burned-trees",my.vector, nl.obj=nl.test1)
NLCommand("show burned-trees", nl.obj=nl.test1)

# submit a list with nested vectors (overrides the variable burned-trees - just for testing, makes no sense)
my.list <- list(c(1,2,3,4,5,6),c(10,11,12),c(TRUE, FALSE, TRUE, FALSE))
print(my.list)
NLCommand("set burned-trees", my.list, nl.obj=nl.test1)
NLCommand("show burned-trees", nl.obj=nl.test1)

# submit a dataframe (overrides the variable burned-trees - just for testing, makes no sense)
my.data.frame <- data.frame(c1=c(1,2,3,4),c2=c(10,11,12,13),c3=c('test1','test2',c4='test3','test4'))
print(my.data.frame)
NLCommand("set burned-trees", my.data.frame, "show burned-trees", nl.obj=nl.test1)


# use NLQuit(nl.obj=nl.test1) to close the NetLogo Window
