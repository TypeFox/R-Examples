library(RNetLogo)

# path to NetLogo installation folder
# PLEASE FILL IN THE PATH TO YOUR NetLogo INSTALLATION FOLDER !!!!
path.to.NetLogo <- "C:/Program Files/NetLogo 5.3/app"
                                                                  
# to start NetLogo without GUI (=headless mode, multiple instances possible)
# and in own, user managed, object instance variable with the name "my.obj1" use: 
my.obj1 <- "my.obj1"
NLStart(path.to.NetLogo, gui=FALSE, nl.obj=my.obj1)
 
# start another NetLogo instance into variable with the name "my.obj2" use: 
my.obj2 <- "my.obj2"
NLStart(path.to.NetLogo, gui=FALSE, nl.obj=my.obj2)

# relative path to a sample model
model.path <- "/models/Sample Models/Earth Science/Fire.nlogo"

# load a sample model into myobj1 instance
NLLoadModel(paste(path.to.NetLogo,model.path,sep=""),nl.obj=my.obj1)

# load a sample model into myobj2 instance
NLLoadModel(paste(path.to.NetLogo,model.path,sep=""),nl.obj=my.obj2)

# set the variable burned-trees with the string "my.obj1 test" of the NetLogo instance my.obj1
NLCommand("set burned-trees \"myobj1 test\"",nl.obj=my.obj1)

# set the variable burned-trees with the string "myobj2 test" of the NetLogo instance my.obj2
NLCommand("set burned-trees \"myobj2 test\"",nl.obj=my.obj2)

name.var1 <- NLReport("burned-trees", nl.obj=my.obj1)
name.var2 <- NLReport("burned-trees", nl.obj=my.obj2)
print(name.var1)
print(name.var2)

# quit both NetLogo instances
NLQuit(my.obj1)
NLQuit(my.obj2)
