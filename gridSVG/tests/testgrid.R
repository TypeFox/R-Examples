# Try to reproduce testsvg.R from grid
library(grid)
library(gridSVG)

dev.new(width=6, height=6)
# Force white background and black foreground
pushViewport(viewport(gp=gpar(col="black", fill="white")))
grid.rect()
# NOTE: that svg.R has transforms that assume that (0, 0)
# is at bottom-left so y-locations and heights in "native"
# coordinates on an X11 device are not
# handled properly (by svg.R;  grid draws them fine).
# So for the outer viewport, we fudge an approximation to the
# "pixel" locations in testsvg.R
# Ultimate fix requires something in svg.R transforms (ty and th)
# Also, the viewport pushed above means that the "native"
# device coordinates are no longer available here
pushViewport(viewport(unit(1, "mm"), unit(1, "mm"),
                       unit(0.5, "npc"), unit(0.6, "npc"),
                       just=c("left", "bottom"),
                       xscale=c(0, 11), yscale=c(0, 11)))
grid.rect(gp=gpar(col="green"))
grid.lines(1:10, 10:1,
           default.units="native",
           gp=gpar(col="green"))
grid.polygon(c(1, 3, 4, 1), c(1, 1, 5, 4),
           default.units="native",
           gp=gpar(fill="grey", col=NA))
grid.rect(rep(6, 2), c(3, 7), 2, 1,
          just=c("left", "bottom"),
          default.units="native",
          gp=gpar(fill="cyan"))
grid.text(c("some text", "some more text!"), 2, 8:7,
          just="left",
          default.units="native")
grid.circle(rep(8, 2), 3, c(.1, 2),
            default.units="native",
            gp=gpar(col="blue", fill=NA))
grid.text("centred text", 4, 5,
          default.units="native", rot=20)
  pushViewport(viewport(x=6, y=5, w=3, h=1,
                         default.units="native",
                         just=c("left", "bottom"),
                         xscale=c(0, 1), yscale=c(0, 1)))
  grid.rect(0, 0, 1, 1,
            just=c("left", "bottom"),
            default.units="native",
            gp=gpar(fill=NA, col="black"))
  grid.text("text in a box", 0.1, 0.5,
            just=c("left", "bottom"),
            default.units="native")
  popViewport()
grid.rect(5, 2, 2, 7,
          default.units="native",
          just=c("left", "bottom"),
          gp=gpar(fill="green", alpha=.5))
popViewport()

popViewport()

grid.export("grid.svg")
dev.off()
