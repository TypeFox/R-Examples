library(grid)
library(gridSVG)

dev.new(width=6, height=6)

grid.newpage()
# Some default settings
pushViewport(viewport(gp=gpar(col="black", fill=NA)))

n <- 30
L <- 50 #mm
A <- 20 #mm
t <- seq(0, 2*pi, length=n)
x <- A*cos(t)
v <- -A*sin(t)
theta <- atan(x/L)
y <- L*cos(theta)

pushViewport(viewport(layout=grid.layout(4, 1,
                         widths=unit(2.2*A, "null"),
                         heights=unit(c(1, L + 10, 2.2*A, 2),
                           c("lines", "null", "null", "lines")),
                         respect=TRUE)))
pushViewport(viewport(layout.pos.row=2,
                       xscale=A * c(-1.1, 1.1),
                       yscale=c(0, L + 10)))
grid.rect()
grid.lines(name="chain",
           x=unit(c(0.5, x[1]), c("npc", "native")),
           y=unit(1, "npc") - unit(c(0, abs(y)[1]), c("npc", "native")))
grid.animate("chain",
             x=animUnit(unit(c(rep(0.5, n), x),
                             c(rep(c("npc", "native"), each=n))),
                        timeid=rep(1:n, 2)),
             y=animUnit(unit.c(unit(rep(1, n), "npc"),
                        unit(1, "npc") - unit(y, "native")),
                        timeid=rep(1:n, 2)),
             duration=5, rep=TRUE)
grid.circle(name="weight",
            x=unit(x[1], "native"),
            y=unit(1, "npc") - unit(y[1], "native"),
            r=unit(1, "mm"),
            gp=gpar(fill="black"))
grid.animate("weight",
             x=unit(x, "native"),
             y=unit(1, "npc") - unit(abs(y), "native"),
             duration=5, rep=TRUE)
popViewport()

pushViewport(viewport(layout.pos.row=3,
                       xscale=A * c(-1.1, 1.1),
                       yscale=A * c(-1.1, 1.1)))
grid.rect()
grid.lines(x=unit(c(0, 0), "native"),
           gp=gpar(lty="dashed", col="grey"))
grid.lines(y=unit(c(0, 0), "native"),
           gp=gpar(lty="dashed", col="grey"))
grid.text("Displacement", y=unit(-1, "lines"))
grid.text("Velocity", x=unit(-1, "lines"), rot=90)
grid.circle(name="key",
            x=unit(x[1], "native"),
            y=unit(v[1], "native"),
            r=unit(1, "mm"),
            gp=gpar(fill="black"))
grid.animate("key",
             x=unit(x, "native"),
             y=unit(v, "native"),
             duration=5, rep=TRUE)
popViewport(2)

popViewport()

grid.export("pendulum.svg")

dev.off()
