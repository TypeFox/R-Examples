
library(grid)
library(gridSVG)

# Primitives to support:
#  rect
grid.newpage()
pushViewport(viewport(width=.5, height=.5, angle=30))
grid.rect()
popViewport()
grid.export("rot-rect-default.svg")
# Check justification
grid.newpage()
pushViewport(viewport(width=.5, height=.5, angle=30))
grid.rect(x=0, y=0, just=c("left", "bottom"))
popViewport()
grid.export("rot-rect-just.svg")

# Primitives to support:
#  text
grid.newpage()
pushViewport(viewport(width=.5, height=.5, angle=30))
grid.text("rot-test")
popViewport()
grid.export("rot-text-default.svg")
# Check justification
grid.newpage()
pushViewport(viewport(width=.5, height=.5, angle=30))
grid.circle(r=unit(1, "mm"))
grid.text("rot-test", just=c("left", "bottom"))
popViewport()
grid.export("rot-text-just.svg")

# Primitives to support:
#  clipPath
grid.newpage()
pushViewport(viewport(width=.5, height=.5, angle=30, clip=TRUE))
grid.circle(r=.6)
popViewport()
grid.export("rot-clip-default.svg")

# Primitives to support:
#  raster
grid.newpage()
pushViewport(viewport(width=.5, height=.5, angle=30))
grid.raster(matrix(1:4/5, ncol=2), interp=FALSE)
popViewport()
grid.export("rot-raster-default.svg")
# Check justification
grid.newpage()
pushViewport(viewport(width=.5, height=.5, angle=30))
grid.raster(matrix(1:4/5, ncol=2), interp=FALSE,
            x=0, y=0, just=c("left", "bottom"))
popViewport()
grid.export("rot-raster-just.svg")

# Primitives to support:
#  plotting symbols
grid.newpage()
pushViewport(viewport(width=.5, height=.5, angle=30))
grid.rect()
grid.points(1:5/6, 1:5/6, pch=1:5)
popViewport()
grid.export("rot-points.svg")

#####################################
# Full STATIC test
library(lattice)
grid.newpage()
pushViewport(viewport(width=.5, height=.5, angle=30))
print(xyplot(mpg ~ disp, mtcars, pch=3),
      newpage=FALSE)
popViewport()
grid.export("rot-lattice.svg")
#####################################

# Implications for animation:
#  rect
grid.newpage()
pushViewport(viewport(width=.5, height=.5, angle=30))
grid.rect(gp=gpar(col="grey"))
grid.rect(x=.1, width=.1, height=.1, just="left", name="r")
grid.animate("r", x=c(.1, .8))
popViewport()
grid.export("rot-rect-animate-x.svg")

grid.newpage()
pushViewport(viewport(width=.5, height=.5, angle=30))
grid.rect(gp=gpar(col="grey"))
grid.rect(x=.1, width=.1, height=.1, just="left", name="r")
grid.animate("r", x=c(.1, .8), height=c(.1, .5))
popViewport()
grid.export("rot-rect-animate-height.svg")

# Implications for animation:
#  text
grid.newpage()
pushViewport(viewport(width=.5, height=.5, angle=30))
grid.rect(gp=gpar(col="grey"))
grid.text("rot-test", x=.1, name="t")
grid.animate("t", x=c(.1, .9))
popViewport()
grid.export("rot-text-animate-x.svg")

# Implications for animation:
#  raster
grid.newpage()
pushViewport(viewport(width=.5, height=.5, angle=30))
grid.rect(gp=gpar(col="grey"))
grid.raster(matrix(1-4:1/20, ncol=2), interpolate=FALSE,
            x=.8, width=.1, just="left")
grid.raster(matrix(1:4/5, ncol=2), interpolate=FALSE,
            x=.1, width=.1, just="left", name="r")
grid.animate("r", x=c(.1, .8))
popViewport()
grid.export("rot-raster-animate-x.svg")

grid.newpage()
pushViewport(viewport(width=.5, height=.5, angle=30))
grid.rect(gp=gpar(col="grey"))
grid.raster(matrix(1-4:1/20, ncol=2), interpolate=FALSE,
            x=.8, width=.1, height=.5, just="left")
grid.raster(matrix(1:4/5, ncol=2), interpolate=FALSE,
            x=.1, width=.1, just="left", name="r")
# NOTE important to specify unit for height because default is NULL!
grid.animate("r", x=c(.1, .8), height=unit(c(.1, .5), "npc"))
popViewport()
grid.export("rot-raster-animate-height.svg")

# Implications for animation:
#  plotting symbols
grid.newpage()
pushViewport(viewport(width=.5, height=.5, angle=30))
grid.rect()
grid.points(5:1/6, 1:5/6, pch=1:5, gp=gpar(col="grey"))
grid.points(1:5/6, 1:5/6, pch=1:5, name="p")
grid.animate("p", x=animUnit(unit(c(1:5, 5:1)/6, "npc"), id=rep(1:5, 2)))
popViewport()
grid.export("rot-points-animate-x.svg")

grid.newpage()
pushViewport(viewport(width=.5, height=.5, angle=30))
grid.rect()
grid.points(5:1/6, 1:5/6, pch=1:5, size=unit(2, "char"), gp=gpar(col="grey"))
grid.points(1:5/6, 1:5/6, pch=1:5, name="p")
grid.animate("p", x=animUnit(unit(c(1:5, 5:1)/6, "npc"), id=rep(1:5, 2)),
             size=animUnit(unit(rep(1:2, each=5), "char"), i=rep(1:5, 2)))
popViewport()
grid.export("rot-points-animate-size.svg")

#####################################
# Full DYNAMIC test
library(lattice)
grid.newpage()
pushViewport(viewport(width=.5, height=.5, angle=30))
print(xyplot(qsec ~ disp, mtcars, pch=3, col="grey",
             ylim=extendrange(c(mtcars$mpg, mtcars$qsec))),
      newpage=FALSE, prefix="plotgrey")
print(xyplot(mpg ~ disp, mtcars, pch=3,
             ylim=extendrange(c(mtcars$mpg, mtcars$qsec))),
      newpage=FALSE, prefix="plot1")
grid.animate("plot1.xyplot.points.panel.1.1", 
             y=animUnit(unit(c(mtcars$mpg, mtcars$qsec), "native"),
                 id=rep(1:nrow(mtcars), 2)))
ylab <- grid.get("plot1.ylab")
grid.animate("plot1.ylab", y=ylab$y - unit(0:1, "in"),
             "fill-opacity"=1:0, "stroke-opacity"=1:0)
grid.animate("plotgrey.ylab", y=ylab$y + unit(1:0, "in"), 
             "fill-opacity"=0:1, "stroke-opacity"=0:1)
popViewport()
grid.export("rot-lattice-animate.svg")
#####################################

# Implications for exported coordinate system info (*.svg.coords.js)
#  need to record rotation info 
#  AND use it in convertViewportX() etc
library(lattice)
grid.newpage()
pushViewport(viewport(width=.5, height=.5, angle=30,
                      name="rot"))
print(xyplot(mpg ~ disp, mtcars, pch=3),
      newpage=FALSE, prefix="coords")
downViewport("coords.panel.1.1.off.vp")
grid.circle(unit(200, "native"), unit(30, "native"),
            r=unit(2, "mm"), gp=gpar(col=NA, fill="grey"))
grid.circle(unit(300, "native"), unit(35, "native"),
            r=unit(1, "mm"), gp=gpar(col=NA, fill="grey"))
upViewport(0)
grid.export("rot-lattice-coords.svg",
            exportCoords="file", exportMappings="file", usePaths="none")
# Read image back into R
library(XML)
svg <- xmlParse("rot-lattice-coords.svg")
# Read in coord info 
gridSVGCoords(readCoordsJS("rot-lattice-coords.svg.coords.js"))
gridSVGMappings(readMappingsJS("rot-lattice-coords.svg.mappings.js"))
# Add new point to panel
panel <- getNodeSet(svg,
                    "//svg:g[contains(@id, 'coords.panel.1.1.off.vp')]",
                    namespaces=c(svg="http://www.w3.org/2000/svg"))[[1]]
vpname <- getSVGMappings("coords.panel.1.1.off.vp", "vp")[1]
pos <- viewportConvertPos(vpname, 200, 30, "native")
circ <- newXMLNode("circle",
                   parent = panel,
                   attrs = list(
                       cx = pos$x, cy = pos$y,
                       r = viewportConvertWidth(vpname, 2, "mm", "svg"),
                       stroke = "red",
                       fill = "red",
                       "fill-opacity" = .5))
dim <- viewportConvertDim(vpname, 100, 5, "native", "svg")
line <- newXMLNode("polyline",
                   parent = panel,
                   attrs = list(
                       points =
                         paste(paste(pos$x, pos$y, sep=","),
                               paste(pos$x + dim$w, pos$y + dim$h, sep=",")),
                       r = viewportConvertWidth(vpname, 2, "mm", "svg"),
                       stroke = "red",
                       fill = "red",
                       "fill-opacity" = .5))
saveXML(svg, "rot-lattice-coords-mod.svg")

# Pre-existing bug in convertViewportWidth() and convertViewportHeight()
library(lattice)
print(xyplot(mpg ~ disp, mtcars, pch=3), prefix="bug")
downViewport("bug.panel.1.1.off.vp")
grid.rect(x=200, y=30, width=100, height=5, just=c("left", "bottom"),
          default="native", gp=gpar(col="grey", fill=NA))
grid.export("bug-lattice-coords.svg",
            exportCoords="file", exportMappings="file", usePaths="none")
library(XML)
svg <- xmlParse("bug-lattice-coords.svg")
gridSVGCoords(readCoordsJS("bug-lattice-coords.svg.coords.js"))
gridSVGMappings(readMappingsJS("bug-lattice-coords.svg.mappings.js"))
panel <- getNodeSet(svg,
                    "//svg:g[contains(@id, 'bug.panel.1.1.off.vp')]",
                    namespaces=c(svg="http://www.w3.org/2000/svg"))[[1]]
vpname <- getSVGMappings("bug.panel.1.1.off.vp", "vp")[1]
x <- viewportConvertX(vpname, 200, "native")
y <- viewportConvertY(vpname, 30, "native")
w <- viewportConvertWidth(vpname, 100, "native", "svg")
h <- viewportConvertHeight(vpname, 5, "native", "svg")
line <- newXMLNode("polyline",
                   parent = panel,
                   attrs = list(points = paste(paste(x, y, sep=","),
                                    paste(x + w, y + h, sep=",")),
                       r = viewportConvertWidth(vpname, 2, "mm", "svg"),
                       stroke = "red",
                       fill = "red",
                       "fill-opacity" = .5))
saveXML(svg, "bug-lattice-coords-mod.svg")

# Implications for exported coordinate system info (*.svg.coords.js)
#  in viewportCreate()
library(lattice)
grid.newpage()
pushViewport(viewport(width=.5, height=.5, angle=30,
                      name="rot"))
print(xyplot(mpg ~ disp, mtcars, pch=3),
      newpage=FALSE, prefix="coords")
downViewport("coords.panel.1.1.off.vp")
grid.yaxis(main=FALSE, gp=gpar(col="grey", fill="grey"))
upViewport(0)
grid.export("rot-lattice-coords-create.svg",
            exportCoords="file", exportMappings="file", usePaths="none")
library(XML)
svg <- xmlParse("rot-lattice-coords-create.svg")
# Read in coord info 
gridSVGCoords(readCoordsJS("rot-lattice-coords-create.svg.coords.js"))
gridSVGMappings(readMappingsJS("rot-lattice-coords-create.svg.mappings.js"))
# Create new viewport
vpname <- getSVGMappings("coords.panel.1.1.off.vp", "vp")[1]
vp <- viewportCreate(vpname)
grid.newpage()
pushViewport(vp)
# Draw an axis and, convert it to SVG, and extract axis SVG content
grid.yaxis(main=FALSE, gp=gpar(col="red", fill="red"))
newsvg <- grid.export(NULL)
axissvg <- getNodeSet(newsvg$svg,
                      "//svg:g[contains(@id, 'yaxis')]",
                      namespaces=c(svg="http://www.w3.org/2000/svg"))[[1]]
# Add new axis to panel
panel <- getNodeSet(svg,
                    "//svg:g[contains(@id, 'coords.panel.1.1.off.vp')]",
                    namespaces=c(svg="http://www.w3.org/2000/svg"))[[1]]
addChildren(panel, kids=list(axissvg))
saveXML(svg, "rot-lattice-coords-create-mod.svg")

# Implications for exported coordinate system info (*.svg.coords.js)
#  in javascript viewportConvertX, etc
library(lattice)
grid.newpage()
pushViewport(viewport(width=.5, height=.5, angle=30,
                      name="rot"))
print(xyplot(mpg ~ disp, mtcars, pch=3),
      newpage=FALSE, prefix="coords")
downViewport("coords.panel.1.1.off.vp")
grid.rect(x=200, y=30, width=100, height=5, just=c("left", "bottom"),
          default="native", gp=gpar(col="grey", fill=NA))
upViewport(0)
grid.garnish("coords.background",
             onclick="addPoint()",
             "pointer-events"="all")
grid.script(file="rot-lattice-coords.js")
grid.export("rot-lattice-coords-js.svg",
            exportCoords="file", exportMappings="file", usePaths="none",
            exportJS="file")

# Pre-existing bug in convertViewportWidth() and convertViewportHeight() (in JS)
library(lattice)
print(xyplot(mpg ~ disp, mtcars, pch=3), prefix="bug")
downViewport("bug.panel.1.1.off.vp")
grid.rect(x=200, y=30, width=100, height=5, just=c("left", "bottom"),
          default="native", gp=gpar(col="grey", fill=NA))
grid.garnish("bug.background",
             onclick="addPoint()",
             "pointer-events"="all")
grid.script(file="bug-lattice-coords.js")
grid.export("bug-lattice-coords-js.svg",
            exportCoords="file", exportMappings="file", usePaths="none",
            exportJS="file")



