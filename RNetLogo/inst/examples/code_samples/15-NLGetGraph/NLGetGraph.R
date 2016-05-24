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
model.path <- "/models/Sample Models/Networks/Preferential Attachment.nlogo"
NLLoadModel(paste(path.to.NetLogo,model.path,sep=""), nl.obj=nl.test1)

# submit a single command
NLCommand("setup", nl.obj=nl.test1)

# run the simulation for four time steps
NLDoCommand(6, "go", nl.obj=nl.test1)

# get the link network as graph
# Package igraph is required! 
my.graph <- NLGetGraph(nl.obj=nl.test1)

# plot the directed network graph
plot(my.graph, layout=layout.kamada.kawai, vertex.label=V(my.graph)$name,
     vertex.shape="rectangle", vertex.size=20, asp=FALSE)

# now, we can execute all gaph/network analysis from igraph package
# for example do shortest path analysis

# calculate the mean shortest path in the network
average.path.length(my.graph, directed=TRUE, unconnected=TRUE)

# get a matrix of shortest paths
shortest.paths(my.graph, v=V(my.graph), mode = c("all", "out", "in"),
      weights = NULL)


# for details on shortest paths see:
?shortest.paths


# one example application:
# test, how the length of mean shortest path increase during simulation:

# reset simulation
NLCommand("setup", nl.obj=nl.test1)

# definition of a function that goes one step ahead and calculates the mean path length
step.ahead <- function(x)
{
  NLCommand("go", nl.obj=nl.test1)
  return (c(x,average.path.length(NLGetGraph(nl.obj=nl.test1), directed=TRUE, unconnected=TRUE)))
}

# call the function "step.aheas" 20 times and save the values of the mean path length in mean.path.length
mean.path.length <- lapply(1:20, function(x) {step.ahead(x)})

# transform list to data.frame
mean.path.length <- data.frame(do.call("rbind",mean.path.length)) 
names(mean.path.length) <- c("time","meanlength")
  
# plot the result
plot(mean.path.length, type='l', main='average path length in network')


# use NLQuit(nl.obj=nl.test1) to close the NetLogo Window


