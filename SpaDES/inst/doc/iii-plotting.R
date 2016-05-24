## ----initial-clearPlot, eval=TRUE, echo=FALSE, message=FALSE-------------
SpaDES::clearPlot()

## ----load_files, eval=TRUE, echo=TRUE, message=FALSE, fig.height=2-------
#  Make list of maps from package database to load, and what functions to use to load them
library(data.table)
library(igraph)
library(raster)
library(SpaDES)
filelist <- data.frame(file = (system.file("maps", package = "SpaDES") %>%
  dir(., full.names = TRUE, pattern = "tif"))[-c(2, 5) ])  # omit forestAge and percentPine for simplicity
print(filelist)
# Load files to memory (using rasterToMemory), assign to a simList we call mySim 
mySim <- loadFiles(filelist = filelist)
# put into a single stack object in the simulation environment for ease of use below
landscape <- stack(mySim$DEM, mySim$forestCover, mySim$habitatQuality)

## ----first_plot, eval=TRUE, echo=TRUE, fig.height=2----------------------
Plot(landscape, new = TRUE)
# make a SpatialPoints object
caribou <- SpatialPoints(coords = cbind(x = stats::runif(1e2, -50, 50),
                                        y = stats::runif(1e2, -50, 50)))
Plot(caribou)
Plot(caribou, addTo = "landscape$habitatQuality")

# from SpatialPolygons help file
Sr1 <- Polygon(cbind(c(2, 4, 4, 1, 2), c(2, 3, 5, 4, 2))*20-50)
Sr2 <- Polygon(cbind(c(5, 4, 2, 5), c(2, 3, 2, 2))*20-50)

Srs1 <- Polygons(list(Sr1), "s1")
Srs2 <- Polygons(list(Sr2), "s2")
SpP <- SpatialPolygons(list(Srs1, Srs2), 1:2)
Plot(SpP)
Plot(SpP, addTo = "landscape$habitatQuality", gp = gpar(lwd = 2))

# from SpatialLines help file
l1 = cbind(c(10, 2, 30), c(30, 2, 2))
l1a = cbind(l1[, 1] + .05, l1[, 2] + .05)
l2 = cbind(c(1, 20, 3), c(10, 1.5, 1))
Sl1 = Line(l1)
Sl1a = Line(l1a)
Sl2 = Line(l2)
S1 = Lines(list(Sl1, Sl1a), ID = "a")
S2 = Lines(list(Sl2), ID = "b")
Sl = SpatialLines(list(S1, S2))
Plot(Sl, gp = gpar(col = c("red", "blue"), lwd = 2), addTo = "landscape$DEM")

## ----mixing_layer_types, eval=TRUE, echo=TRUE, fig.height=5--------------
Plot(landscape, caribou, mySim$DEM, SpP, new = TRUE, axes = TRUE,
     gp = gpar(cex = 0.5), visualSqueeze = 0.7)

## ----ggplot, eval=TRUE, echo=TRUE, cache=TRUE, fig.height=2--------------
library(ggplot2)
ggObj <- qplot(stats::rnorm(1e3), binwidth = 0.1)

Plot(caribou, ggObj, new = TRUE) 

## ----simList, eval=TRUE, echo=TRUE, cache=TRUE, fig.height=2-------------
Plot(mySim, new = TRUE) 

## ----set_colours, eval=TRUE, echo=TRUE, fig.height=2---------------------
library(RColorBrewer)

# can change colour palette
Plot(landscape, new = TRUE) # original

mapColours <- list(
  DEM = topo.colors(50),
  forestCover = colorRampPalette(c("blue", "orange", "purple", "red"))(50),
  habitatQuality = brewer.pal(9, "Spectral")
)
setColors(landscape, n = 50) <- mapColours
Plot(landscape, new = TRUE) # oh, how pretty!

## ----gp_gpAxis_gpText, eval=TRUE, echo=TRUE, fig.height=2----------------
Plot(caribou, new = TRUE, gpAxis = gpar(cex = 0.4), size = 1)
Plot(mySim$DEM, gpText = gpar(cex = 0.4))
Plot(mySim$DEM, caribou, gpText = list(gpar(cex = 2), gpar(cex = 0.1)), new = TRUE)

## ----visualSqueeze, eval=TRUE, echo=TRUE, fig.height=2-------------------
# x axis gets cut off in pdf and html
Plot(mySim$DEM, new = TRUE)
Plot(mySim$DEM, visualSqueeze = 0.6, new = TRUE)

## ----legendRange, eval=TRUE, echo=TRUE, fig.height=2---------------------
Plot(mySim$DEM, legendRange = c(0, 500), new = TRUE)

## ----zoomExtent, eval=TRUE, echo=TRUE, fig.height=2----------------------
Plot(mySim$DEM, zoomExtent = extent(c(-1, 10, -1, 20)), new = TRUE)

## ----arrows, eval=TRUE, echo=TRUE, fig.height=2--------------------------
Plot(mySim$DEM, new = TRUE)
Plot(Sl, addTo = "mySim$DEM", length = 0.2)

## ----simple_add, eval=TRUE, echo=TRUE, fig.height=3----------------------
Plot(landscape, new = TRUE)
# can add a new plot to the plotting window
Plot(caribou, new = FALSE, axes = FALSE)

## ----add_with_rearrangement, eval=TRUE, echo=TRUE, fig.height=2----------
Plot(landscape, new = TRUE)
# can add a new plot to the plotting window
Plot(caribou, new = FALSE, axes = FALSE)

## ----add_with_same_name, eval=TRUE, echo=TRUE, fig.height=2--------------
Plot(landscape, new = TRUE)
landscape$forestCover[] = ((landscape$forestCover[] +10) %% 30)
# can add a new plot to the plotting window
Plot(landscape, new = FALSE)
# note: zeros are treated as no colour by default.
# if this is not the correct behavior, use `zero.color=NULL`
Plot(landscape, new = FALSE, zero.color = NULL)

## ----speedup, eval=TRUE, echo=TRUE, fig.height=2-------------------------
system.time(Plot(landscape, caribou, mySim$DEM, new = TRUE))
system.time(Plot(landscape, caribou, mySim$DEM, speedup = 10, new = TRUE))
# can add a new plot to the plotting window

## ----add, eval=TRUE, echo=TRUE, fig.height=2-----------------------------
Plot(landscape, new = TRUE)
Plot(caribou, addTo = "landscape$DEM", size = 2, axes = FALSE)

## ----clearPlot, eval=TRUE, echo=TRUE, fig.height=2-----------------------
clearPlot()
Plot(caribou)

## ----clickValues, eval=FALSE, echo=TRUE----------------------------------
#  Plot(landscape, new = TRUE)
#  clickValues(3) # click at three locations on the Plot device

## ----clickExtent, eval=FALSE, echo=TRUE----------------------------------
#  Plot(landscape, new = TRUE)
#  clickExtent() # click at two locations on the Plot device

## ----rePlot, eval=FALSE, echo=TRUE, cache=TRUE---------------------------
#  rePlot()
#  rePlot(4)
#  rePlot(visualSqueeze = 1, axes = FALSE)

## ----Plot a .spadesPlot object, eval=FALSE, echo=TRUE, cache=TRUE--------
#  plots <- Plot(landscape, new = TRUE)
#  
#  # change values
#  landscape$forestCover[landscape$habitatQuality > 0.9] <- 0
#  
#  Plot(plots)
#  # same as:
#  rePlot()
#  
#  # but can be combined with other objects
#  Plot(caribou, plots, new = TRUE)
#  Plot(caribou, plots, Sl)

