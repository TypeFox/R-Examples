## ----install-fastshp, echo=TRUE, eval=FALSE------------------------------
#  install.packages("fastshp", repos="http://rforge.net", type="source")

## ----SpaDES-demo, eval=FALSE, echo=TRUE----------------------------------
#  library("SpaDES")
#  demo("spades-simulation", package="SpaDES")

## ----SpaDES-modules, eval=FALSE, echo=TRUE-------------------------------
#  saveTo <- "~/SpaDES-modules" # change this to suit your needs
#  downloadModule(name="moduleName", path=saveTo)

## ----view-sim, eval=FALSE, echo=TRUE-------------------------------------
#  # full simulation details:
#  #  simList object info + simulation data
#  mySim
#  
#  # less detail:
#  # simList object isn't shown; object details are
#  ls.str(mySim)
#  
#  # least detail:
#  # simList object isn't shown; object names only
#  ls(mySim)

## ----view-dependencies, eval=FALSE, echo=TRUE----------------------------
#  library(igraph)
#  depsEdgeList(mySim, FALSE)  # data.frame of all object dependencies
#  moduleDiagram(mySim)        # plots simplified module (object) dependency graph
#  objectDiagram(mySim)        # plots object dependency diagram

## ----view-event-sequences, eval=FALSE, echo=TRUE-------------------------
#  options(spades.nCompleted = 50)   # default: store 10 events in the completed event list
#  mySim <- simInit(...)             # initialize a simulation using valid parameters
#  mySim <- spades(mySim)            # run the simulation, returning the completed sim object
#  eventDiagram(mySim, "0000-06-01") # visualize the sequence of events for all modules

## ----download-module, echo=TRUE, eval=FALSE------------------------------
#  downloadModule("moduleName", "localpath/to/put/my/modules/directory")

