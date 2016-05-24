library(grid)
library(gridSVG)

# animValues

animValue(letters[1:4])
animValue(letters[1:4], timeid=rep(1:2, 2))
animValue(letters[1:4], id=rep(1:2, 2))

as.animValue(letters[1:4])
as.animValue(matrix(letters[1:12], ncol=4))
as.animValue(matrix(letters[1:12], ncol=4), multVal=TRUE)
as.animValue(list(letters[1:3], letters[4:6]))
as.animValue(list(letters[1:3], letters[4:6]), multVal=TRUE)

# animUnits

animUnit(unit(1:4, "npc"))
animUnit(unit(1:4, "npc"), timeid=rep(1:2, 2))
animUnit(unit(1:4, "npc"), id=rep(1:2, 2))

as.animUnit(1:4, "npc")
as.animUnit(unit(1:4, "npc"))
as.animUnit(matrix(1:12, ncol=4), "in")
as.animUnit(matrix(1:12, ncol=4), "in", multVal=TRUE)
as.animUnit(list(unit(1:3, "npc"), unit(4:6, "in")))
as.animUnit(list(unit(1:3, "npc"), unit(4:6, "in")), multVal=TRUE)

# Some default settings
pushViewport(viewport(gp=gpar(col="black", fill=NA)))

grid.rect(name="rect",
          x=0,
          y=0,
          just=c("left", "bottom"))
grid.animate("rect", x=unit(0:30, "mm"), duration=5, rep=TRUE)
grid.circle(name="circle",
            x=unit(0.5, "npc") + unit(0, "mm"),
            r=unit(10, "mm"))
grid.animate("circle", x=unit(0.5, "npc") + unit(0:30, "mm"),
             duration=5, rep=TRUE)
grid.text("hello", name="text1",
          x=unit(0.3, "npc") + unit(0, "mm"))
grid.animate("text1",
             x=unit(0.3, "npc") + unit(0:30, "mm"),
             duration=5, rep=TRUE)
grid.text("hello", name="text2",
          x=unit(0.3, "npc") + unit(0, "mm"),
          y=unit(0.3, "npc") + unit(0, "mm"))
grid.animate("text2",
             x=unit(0.3, "npc") + unit(0:30, "mm"),
             y=unit(0.3, "npc") + unit(0:30, "mm"),
             duration=5, rep=TRUE)

popViewport()

grid.export("animate.svg")


# Animating rectangles

# There are numerous possibilities to consider:
#   The animation values could be numeric, unit, matrix, or list
#   The original values could spec a single rect or multiple rects
#   We could animate only one of x/y/width/height or several of them at once

# Simple case
# (single rect, anim only x, anim values are just numeric)
grid.newpage()
grid.text("One rectangle moves across",
          y=unit(1, "lines"))
grid.rect()
grid.rect(x=.2, y=.2, width=.1, height=.1, name="rect")
grid.animate("rect", x=c(.2, .8), duration=3)
grid.export("anim-rect-simple.svg")

# Complex case
# (multiple rects, anim x/y/width/height, anim values are matrices and lists)
grid.newpage()
grid.text("Three rectangles: one goes up, one goes across, and
one goes diagonal and gets smaller",
          y=unit(1, "lines"))
grid.rect()
grid.rect(x=rep(.2, 3), y=.2, width=.1, height=.1, name="rect")
grid.animate("rect",
             x=cbind(c(.2, .8), c(.2, .8), .2),
             y=cbind(.2, c(.2, .8), c(.2, .8)),
             width=list(unit(.1, "npc"),
               unit(c(.1, 1), c("npc", "cm")),
               unit(.1, "npc")),
             height=list(unit(.1, "npc"),
               unit(c(.1, 1), c("npc", "cm")),
               unit(.1, "npc")),             
             duration=3)
grid.export("anim-rect-complex.svg")

# Animating circles

# Complex case
# (multiple circles, anim x/y/width/height, anim values are matrices and lists)
grid.newpage()
grid.text("Three circles: one goes up, one goes across, and
one goes diagonal and gets smaller",
          y=unit(1, "lines"))
grid.rect()
grid.circle(x=rep(.2, 3), y=.2, r=.1, name="circle")
grid.animate("circle",
             x=cbind(c(.2, .8), c(.2, .8), .2),
             y=cbind(.2, c(.2, .8), c(.2, .8)),
             r=list(unit(.1, "npc"),
               unit(c(.1, 1), c("npc", "cm")),
               unit(.1, "npc")),             
             duration=3)
grid.export("anim-circle-complex.svg")

# Animating points

# Complex case
# (multiple circles, anim x/y/width/height, anim values are matrices and lists)
grid.newpage()
grid.text("Three points: one goes up, one goes across, and
one goes diagonal and gets larger",
          y=unit(1, "lines"))
grid.rect()
pushViewport(viewport())
grid.points(x=rep(.2, 3), y=rep(.2, 3), size=unit(2, "mm"), name="points")
grid.animate("points",
             x=cbind(c(.2, .8), c(.2, .8), .2),
             y=cbind(.2, c(.2, .8), c(.2, .8)),
             size=list(unit(2, "mm"),
               unit(c(2, .1), c("mm", "npc")),
               unit(2, "mm")),             
             duration=3)
grid.export("anim-points-complex.svg")

# Animating text

# Complex case
# (multiple text, anim x/y/width/height, anim values are matrices and lists)
grid.newpage()
grid.text("Three letters:  one goes up, one goes across, and
one goes diagonal",
          y=unit(1, "lines"))
grid.rect()
grid.text(letters[1:3], x=rep(.2, 3), y=.2, name="text")
grid.animate("text",
             x=cbind(c(.2, .8), c(.2, .8), .2),
             y=cbind(.2, c(.2, .8), c(.2, .8)),             
             duration=3)
grid.export("anim-text-complex.svg")

# Animating lines

# Simple case
# (line only has two points, animation only has two points, only animate x)
grid.newpage()
grid.text("45 degree line becomes vertical",
          y=unit(1, "lines"))
grid.rect()
grid.lines(c(.1, .9), c(.1, .9), name="lines")
grid.animate("lines",
             x=cbind(c(.1, .9), c(.5, .5)),
             duration=3)
grid.export("anim-lines-simple.svg")

# Complex case
# (line has many points, animation has three points, only animate y)
x <- seq(-pi, pi, length.out=100)
y <- sin(x)
grid.newpage()
grid.text("Sine curve becomes flat then inverts (on y)",
          y=unit(1, "lines"))
grid.rect()
pushViewport(dataViewport(x, y))
grid.lines(x, y, default.units="native", name="lines")
grid.animate("lines",
             y=cbind(y, 0, -y),
             duration=3)
grid.export("anim-lines-complex.svg")

# Animating polylines

# Simple case
# (line only has two points, animation only has two points, only animate x)
grid.newpage()
grid.text("Two parallel lines slide to the right",
          y=unit(1, "lines"))
grid.rect()
grid.polyline(c(.1, .2, .3, .4),
              c(.1, .9, .1, .9),
              id=rep(1:2, each=2), name="polyline")
grid.animate("polyline",
             x=animUnit(unit(c(.1, .2, .3, .4,
                               .5, .6, .7, .8),
                             unit="npc"),
                           id=rep(rep(1:2, each=2), 2),
                           timeid=rep(1:2, each=4)),
             duration=3)
grid.export("anim-polyline-simple.svg")

# Complex case
# (line only has many points, animation only has many points, animate x and y)
grid.newpage()
grid.text("Two random walks",
          y=unit(1, "lines"))
grid.rect()
n <- 50
x <- 1:n
set.seed(1000)
y1 <- runif(n, .6, .8)
y2 <- runif(n, .2, .4)
pushViewport(dataViewport(x, yscale=0:1))
grid.polyline(rep(x[1:2], 2), c(y1[1:2], y2[1:2]),
              default.units="native",
              id=rep(1:2, each=2), name="polyline")
grid.animate("polyline",
             x=animUnit(unit(rep(x[unlist(lapply(2:n, seq))], 2),
                             "native"),
                        id=rep(1:2, each=sum(2:n)),
                        timeid=rep(1:(n - 1), 2:n)),
             y=animUnit(unit(c(y1[unlist(lapply(2:n, seq))],
                               y2[unlist(lapply(2:n, seq))]),
                             "native"),
                        id=rep(1:2, each=sum(2:n)),
                        timeid=rep(1:(n - 1), 2:n)),
             duration=10)
grid.export("anim-polyline-complex.svg")

# Animating segments

# Simple case
# (single segment, animation only has two values, only animate x0)
grid.newpage()
grid.text("45 degree line becomes vertical (on right)",
          y=unit(1, "lines"))
grid.rect()
grid.segments(.1, .1, .9, .9, name="segments")
grid.animate("segments",
             x0=c(.1, .9),
             duration=3)
grid.export("anim-segments-simple.svg")

# Complex case
# (multiple segments, animation has three values, animate x0 and y0)
grid.newpage()
grid.text("crossed lines swing out to vertical then shorten",
          y=unit(1, "lines"))
grid.rect()
grid.segments(c(.1, .9), .1, c(.9, .1), .9, name="segments")
grid.animate("segments",
             x0=cbind(c(.1, .9, .9), c(.9, .1, .1)),
             y0=c(.1, .1, .5),
             duration=3)
grid.export("anim-segments-complex.svg")

# Animating polygons

# Simple case
# (polygon only has three points,
#  animation only has two points, only animate x)
grid.newpage()
grid.text("Single polygon slides to the right",
          y=unit(1, "lines"))
grid.rect()
grid.polygon(c(.1, .2, .3),
              c(.4, .6, .4), name="polygon")
grid.animate("polygon",
             x=animUnit(unit(c(.1, .2, .3,
                               .7, .8, .9),
                             unit="npc"),
                        timeid=rep(1:2, each=3)),
             duration=3)
grid.export("anim-polygon-simple.svg")

# Complex case
# (two polygons, animation has many points, animate x and y)
grid.newpage()
grid.text("Two polygons shrink and grow (flipped) then revert",
          y=unit(1, "lines"))
grid.rect()
grid.polygon(c(.2, .3, .4,
               .6, .7, .8),
             c(.4, .6, .4, .6, .4, .6),
             id=rep(1:2, each=3), name="polygon")
grid.animate("polygon",
             x=animUnit(unit(c(.2, .3, .4,
                               .4, .3, .2,
                               .2, .3, .4,
 
                               .6, .7, .8,
                               .8, .7, .6,
                               .6, .7, .8),
                             "npc"),
                        id=rep(1:2, each=9),
                        timeid=rep(rep(1:3, each=3), 2)),
             y=animUnit(unit(c(.4, .6, .4,
                               .6, .4, .6,
                               .4, .6, .4,
                               
                               .6, .4, .6,
                               .4, .6, .4,
                               .6, .4, .6),
                             "npc"),
                        id=rep(1:2, each=9),
                        timeid=rep(rep(1:3, each=3), 2)),
             duration=5)
grid.export("anim-polygon-complex.svg")

# Animating paths

# Simple case
# (path has one sub-path,
#  animation only has two points, only animate x)
grid.newpage()
grid.text("Single simple path (triangle) slides to the right",
          y=unit(1, "lines"))
grid.rect()
grid.path(c(.1, .2, .3),
          c(.4, .6, .4),
          gp=gpar(fill="black"),
          name="path")
grid.animate("path",
             x=animUnit(unit(c(.1, .2, .3,
                               .7, .8, .9),
                             unit="npc"),
                        timeid=rep(1:2, each=3)),
             duration=3)
grid.export("anim-path-simple.svg")

# Complex case
# (two polygons, animation has many points, animate x and y)
grid.newpage()
grid.text("Single complex path transmogrifies as it slides to the right",
          y=unit(1, "lines"))
grid.rect()
grid.path(c(.1, .1, .4, .4,
            .2, .2, .3, .3),
          c(.2, .8, .8, .2,
            .4, .6, .6, .4),
          id=rep(1:2, each=4),
          rule="evenodd",
          gp=gpar(fill="black"),
          name="path")
grid.animate("path",
             x=animUnit(unit(c(.1, .1, .4, .4,
                               .2, .2, .3, .3,

                               .35, .35, .65, .65,
                               .45, .45, .55, .55,
               
                               .6, .6, .9, .9,
                               .7, .7, .8, .8),
                             unit="npc"),
                        id=rep(rep(1:2, each=4), 3),
                        timeid=rep(1:3, each=8)),
             y=animUnit(unit(c(.2, .8, .8, .2,
                               .4, .6, .6, .4,

                               .4, .6, .6, .4,
                               .2, .8, .8, .2,

                               .2, .8, .8, .2,
                               .4, .6, .6, .4),
                             unit="npc"),
                        id=rep(rep(1:2, each=4), 3),
                        timeid=rep(1:3, each=8)),               
             duration=3)
grid.export("anim-path-complex.svg")

# Simple case
# (single raster, anim only x, anim values are just numeric)
grid.newpage()
grid.text("One raster moves across",
          y=unit(1, "lines"))
grid.rect()
grid.raster(1:10/11, x=.2, y=.2, width=.1, height=.1, name="raster")
grid.animate("raster", x=c(.2, .8), duration=3)
grid.export("anim-raster-simple.svg")

# Complex case
# (multiple rasters, anim x/y/width/height, anim values are matrices and lists)
grid.newpage()
grid.text("Three rasters: one goes up, one goes across, and
one goes diagonal and gets smaller",
          y=unit(1, "lines"))
grid.rect()
grid.raster(1:10/11, x=rep(.2, 3), y=.2, width=.1, height=.1, name="raster")
grid.animate("raster",
             x=cbind(c(.2, .8), c(.2, .8), .2),
             y=cbind(.2, c(.2, .8), c(.2, .8)),
             width=list(unit(.1, "npc"),
               unit(c(.1, 1), c("npc", "cm")),
               unit(.1, "npc")),
             height=list(unit(.1, "npc"),
               unit(c(.1, 1), c("npc", "cm")),
               unit(.1, "npc")),             
             duration=3)
grid.export("anim-raster-complex.svg")

# Simple case
# (single xspline, anim only x, anim values are just numeric)
grid.newpage()
grid.text("Two xsplines move across",
          y=unit(1, "lines"))
grid.rect()
grid.xspline(c(.3, .1, .5, .3),
             c(.2, .5, .5, .2),
             open=TRUE, shape=1,
             name="xspline-1")
grid.xspline(c(.3, .1, .5),
             c(.6, .9, .9),
             open=FALSE, shape=1,
             gp=gpar(fill="grey"),
             name="xspline-2")
grid.animate("xspline-1",
             x=animUnit(unit(c(.3, .1, .5, .3,
                               .7, .5, .9, .7),
                             "npc"),
                        timeid=rep(1:2, each=4)),
             duration=3)
grid.animate("xspline-2",
             x=cbind(c(.3, .1, .5),
                     c(.7, .5, .9)),
             duration=3)
grid.export("anim-xspline-simple.svg")

# Complex case
# (four xsplines, animation has many points, animate x and y)
grid.newpage()
grid.text("Four xsplines shrink and grow (flipped) then revert",
          y=unit(1, "lines"))
grid.rect()
grid.xspline(c(.3, .1, .5, .3,
               .7, .5, .9, .7),
             c(.2, .5, .5, .2,
               .2, .5, .5, .2),
             shape=1,
             id=rep(1:2, each=4), name="xspline-open")
grid.xspline(c(.3, .1, .5,
               .7, .5, .9),
             c(.6, .9, .9,
               .6, .9, .9),
             open=FALSE, shape=1,
             gp=gpar(fill="grey"),
             id=rep(1:2, each=3), name="xspline-closed")
grid.animate("xspline-open",
             x=animUnit(unit(c(.3, .1, .5, .3,
                               .5, .3, .7, .5,
                               .7, .5, .9, .7,

                               .7, .5, .9, .7,
                               .5, .3, .7, .5,
                               .3, .1, .5, .3),
                             "npc"),
                        id=rep(1:2, each=12),
                        timeid=rep(rep(1:3, each=4), 2)),
             y=animUnit(unit(c(.2, .5, .5, .2,
                               .5, .2, .2, .5,
                               .2, .5, .5, .2,
               
                               .2, .5, .5, .2,
                               .5, .2, .2, .5,
                               .2, .5, .5, .2),
                             "npc"),
                        id=rep(1:2, each=12),
                        timeid=rep(rep(1:3, each=4), 2)),
             duration=5)
grid.animate("xspline-closed",
             x=animUnit(unit(c(.3, .1, .5,
                               .5, .3, .7,
                               .7, .5, .9,

                               .7, .5, .9,
                               .5, .3, .7,
                               .3, .1, .5),
                             "npc"),
                        id=rep(1:2, each=9),
                        timeid=rep(rep(1:3, each=3), 2)),
             y=animUnit(unit(c(.6, .9, .9,
                               .9, .6, .6,
                               .6, .9, .9,

                               .6, .9, .9,
                               .9, .6, .6,
                               .6, .9, .9),
                             "npc"),
                        id=rep(1:2, each=9),
                        timeid=rep(rep(1:3, each=3), 2)),
             duration=5)
grid.export("anim-xspline-complex.svg")


############################################

# Multiple animations on same grob
grid.newpage()
grid.rect(x=.1, y=.1, width=.1, height=.1, name="r")
grid.animate("r", x=c(.1, .9))
grid.animate("r", x=c(.9, .1), begin=3)
grid.export("anim-rect-multi.svg")
 
# Animate group
grid.newpage()
grid.rect(x=.1, y=.1, width=.1, height=.1, name="r")
grid.animate("r", visibility=c("visible", "hidden"), group=TRUE)
grid.export("anim-group.svg")

