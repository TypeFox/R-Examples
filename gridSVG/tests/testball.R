library(grid)
library(lattice)
library(gridSVG)

postscript(width=8, height=6, paper="special")

# Some default settings
pushViewport(viewport(gp=gpar(col="black", fill=NA)))

y <- 1:4
x <- 1:4
g <- factor(c("Earth", "Moon", "Jupiter", "Mars"))

n <- 30
d <- 20 # metres
k <- c(9.8, 1.6, 24.8, 3.7) # gravities
times <- 2 * sqrt(d/k) # (twice) the time taken to fall d metres

cx <- unit(rep(0.5, 2*n - 1), "npc")
dy <- 10*seq(0, 2, length=n)^2
ecy <- unit(40 - c(dy, rev(dy)[-1]), "native")
  
ballpanel <- function(x, y, subscripts) {
  pushViewport(viewport(yscale=c(-10, 50)))
  grid.rect(y=unit(0, "npc"), 
            height=unit(10, "native"),
            just="bottom",
            gp=gpar(fill="grey"))
  duration <- switch(subscripts,
                     times[1],
                     times[2],
                     times[3],
                     times[4])
  col <- switch(subscripts,
                "blue", "grey", "brown", "red")
  grid.circle(name=col,
              x=cx[1],
              y=ecy[1],
              # r=unit(2, "mm"),
              r=unit(1, "native"),
              gp=gpar(col="black", fill=col))
  grid.animate(col, x=cx, y=ecy, duration=duration, rep=TRUE)
  if (subscripts == 1) {
    grid.text("20 metres", 
              x=unit(-1, "lines"), 
              y=unit(20, "native"), 
              just="bottom", rot=270)
    grid.lines(x=unit(-1, "lines"),
               y=unit.c(unit(0, "native"),
                        unit(20, "native") - 
                        unit(0.5, "strwidth", "40 metres") - 
                        unit(2, "mm")))
    grid.lines(x=unit(c(-.75, -1.25), "lines"),
               y=unit(0, "native"))
    grid.lines(x=unit(-1, "lines"),
               y=unit.c(unit(40, "native"),
                        unit(20, "native") + 
                        unit(0.5, "strwidth", "40 metres") + 
                        unit(2, "mm")))
    grid.lines(x=unit(c(-.75, -1.25), "lines"),
               y=unit(40, "native"))
  }
  popViewport()
}

print(xyplot(y ~ x | g, subscripts=TRUE,
             layout=c(4, 1),
             xlab=NULL, ylab=NULL,
             panel=ballpanel,
             scales=list(draw=FALSE)),
      position=c(0.1, 0.1, 0.9, 0.9),
      newpage=FALSE)

popViewport()

grid.export("ball.svg")

dev.off()
