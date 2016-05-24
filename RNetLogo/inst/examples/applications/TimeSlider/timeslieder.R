
# simple example without NetLogo
library(rpanel)

timedata.list1 <- list(rnorm(100), rnorm(100), rnorm(100), rnorm(100))

plottimedata <- function(timedata.list,...) {
   timeslider.draw <- function(panel) {
     plot(timedata.list[[panel$t]], ...)
     mtext(paste("at time ",panel$t),3)
     panel
     }
   panel <- rp.control()
   rp.slider(panel, resolution=1, var=t, from=1, to=length(timedata.list), title="Time", showvalue=TRUE, action = timeslider.draw) }

plottimedata(timedata.list=timedata.list1, main="My title", xlab="bla", ylab="blub")


####################################################################
# usage example with NetLogo
####################################################################
library(RNetLogo)


plottimedata <- function(timedata.list,x.var,y.var,color.var=NULL, ...) {
   timeslider.draw <- function(panel) {
     index.x <- which(names(timedata.list[[panel$t]])==x.var)
     index.y <- which(names(timedata.list[[panel$t]])==y.var)
     color <- "black"
     if (!is.null(color.var)) {
      index.color <- which(names(timedata.list[[panel$t]])==color.var)
      color <- timedata.list[[panel$t]][[index.color]]
      color[color==F] <- "green"
      color[color==T] <- "red"
     }
     plot(timedata.list[[panel$t]][[index.x]],timedata.list[[panel$t]][[index.y]], col=color, ...)
     mtext(paste("at time ",panel$t),3)
     panel
     }
   panel <- rp.control()
   rp.slider(panel, resolution=1, var=t, from=1, to=length(timedata.list), title="Time", showvalue=TRUE, action = timeslider.draw)
}


nl.path <- "C:/Program Files/NetLogo 5.3/app"
model.path <- "/models/Sample Models/Biology/Tumor.nlogo"
NLStart(nl.path)
NLLoadModel(paste(nl.path,model.path,sep=""))
NLCommand("setup")

nruns <- 40
timedata <- list()

for(i in 1:nruns) {
  NLCommand("go")
  timedata[[i]] <- NLGetAgentSet(c("xcor","ycor","who","stem?","metastatic?"), "turtles")
}

world.dim <- NLReport(c("(list min-pxcor max-pxcor)", "(list min-pycor max-pycor)"))
colors <- c("green","red")

plottimedata(timedata.list=timedata, x.var="xcor", y.var="ycor", xlab="x", ylab="y", color.var="metastatic?", main="Tumor cells", xlim=world.dim[[1]] ,ylim=world.dim[[2]])

# close NetLogo
NLQuit()