library(RNetLogo)
library(rpanel)

# define a function to set a logical variable to colors
color.func <- function(color.var,colors,timedata) {
  color <- NULL
  if (!is.null(color.var)) {
    index.color <- which(names(timedata)==color.var)
    color <- timedata[[index.color]]
    color[color==F] <- colors[1]
    color[color==T] <- colors[2]
  }
  return(color)
}

# define a custom plot function using rp.slider of the rpanel package
# to browse through the plots
plottimedata <- function(timedata.list,x.var,y.var,boxplot.var1,
                         boxplot.var2,color.var1=NULL,colors1="black",
                         color.var2=NULL,colors2="black", mains = NULL, ...) {
   # the drawing function, called when the slider position is changed
   timeslider.draw <- function(panel) {
     index.x <- which(names(timedata.list[[panel$t]])==x.var)
     index.y <- which(names(timedata.list[[panel$t]])==y.var)
     index.b1 <- which(names(timedata.list[[panel$t]])==boxplot.var1)
     index.b2 <- which(names(timedata.list[[panel$t]])==boxplot.var2)

     # if a color variable (logical) is given set the colors
     color1 <- color.func(color.var1,colors1,timedata.list[[panel$t]])
     color2 <- color.func(color.var2,colors2,timedata.list[[panel$t]])

     # 4 figures arranged in 2 rows and 2 columns with one title text line
     par(mfrow=c(2,2),oma = c( 0, 0, 1, 0 ))
     # create current plot
     plot(timedata.list[[panel$t]][[index.x]],
          timedata.list[[panel$t]][[index.y]], col=color1, main=mains[1], ...)
     plot(timedata.list[[panel$t]][[index.x]],
          timedata.list[[panel$t]][[index.y]], col=color2, main=mains[2], ...)
     boxplot(timedata.list[[panel$t]][[index.b1]], main=mains[3])
     boxplot(timedata.list[[panel$t]][[index.b2]], main=mains[4])
     title( paste("at time ",panel$t), outer = TRUE )
     panel
   }
   # create a control panel (hosting the slider)
   panel <- rp.control()
   # create a slider to switch the plot data
   rp.slider(panel, resolution=1, var=t, from=1, to=length(timedata.list),
             title="Time", showvalue=TRUE, action = timeslider.draw)
}

# initialize NetLogo
nl.path <- "C:/Program Files/NetLogo 5.3/app"
model.path <- "/models/Sample Models/Biology/Virus.nlogo"
NLStart(nl.path)
# load the Tumor model
NLLoadModel(paste(nl.path,model.path,sep=""))
# initialize the model
NLCommand("setup")
# run the model for 100 time steps and save the turtles of
# every step in one entry of the timedata list
nruns <- 100
timedata <- list()
for(i in 1:nruns) {
  NLCommand("go")
  timedata[[i]] <- NLGetAgentSet(c("who","xcor","ycor","age","sick?","immune?","sick-count"),
                                 "turtles")
}
# get the world dimension to use for the plot
world.dim <- NLReport(c("(list min-pxcor max-pxcor)",
                        "(list min-pycor max-pycor)"))
# define colors to be used for turtle visualization
colors1 <- c("green","red")
colors2 <- c("red","green")
# call the plottimedata function to brwose through the timedata list
plottimedata(timedata.list=timedata, x.var="xcor", y.var="ycor", xlab="x",
             ylab="y", color.var1="sick?", color.var2="immune?",
             boxplot.var1="sick-count", boxplot.var2="age",
             colors1=colors1, colors2=colors2,
             mains=c("Sick","Immune","Stick-count","Age"),
             xlim=world.dim[[1]], ylim=world.dim[[2]])
             
# close NetLogo     
NLQuit()