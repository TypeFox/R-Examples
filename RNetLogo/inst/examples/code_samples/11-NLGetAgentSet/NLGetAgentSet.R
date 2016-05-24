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
path.to.local.file <- 'C:/Users/jthiele/Documents/R/win-library/3.2/RNetLogo/examples/code_samples/11-NLGetAgentSet/agentsettest.nlogo'
NLLoadModel(path.to.local.file, nl.obj=nl.test1) 

# submit a single command
NLCommand("setup", nl.obj=nl.test1)

# NLGetAgentSet is very flexible: patches, turtles and links can be requested. In case of patches it's equivalent to NLGetPatches!
# But in general, everything is also possible with complex requests via NLReport!

# get informations of all breeds "myturtles1", always sorted by "who" variable
my.turtles.1 <- NLGetAgentSet(c("who","xcor","ycor","color"), "myturtles1", nl.obj=nl.test1)
print(my.turtles.1)

# same as before but as list (where each list element represent a agent variable)
my.turtles.1.list <- NLGetAgentSet(c("who","xcor","ycor","color"), "myturtles1", as.data.frame=FALSE, nl.obj=nl.test1)
print(my.turtles.1.list)

# same as before but as a list where each list element represent an agent 
# (ATTENTION: this is very slow! Use it carefully, especially on a large number of agents.)
my.turtles.1.list2 <- NLGetAgentSet(c("who","xcor","ycor","color"), "myturtles1", as.data.frame=FALSE, agents.by.row=TRUE, nl.obj=nl.test1)
print(my.turtles.1.list2)

# get informations from a subset of "myturtles1" with xcor < 3
subset.my.turtles.1 <- NLGetAgentSet(c("who","xcor","ycor","color"), "myturtles1 with [xcor < 3]", nl.obj=nl.test1)
print(subset.my.turtles.1)

# get "who" of the ends of links
links <- NLGetAgentSet(c("[who] of end1","[who] of end2"), "links", nl.obj=nl.test1)
print(links)

# get all patches which has a exactly one of breed "myturtles2" on it, always sorted from upper left to lower right
patchset <- NLGetAgentSet(c("pxcor","pycor","plabel","one-of [who] of myturtles2-here"), "patches with [count myturtles2-here = 1]", nl.obj=nl.test1)
print(patchset)

# use NLQuit(nl.obj=nl.test1) to close the NetLogo Window
