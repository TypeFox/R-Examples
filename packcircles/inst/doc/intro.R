## ------------------------------------------------------------------------

ncircles <- 200
limits <- c(-50, 50)
inset <- diff(limits) / 3
rmax <- 20

xyr <- data.frame(
  x = runif(ncircles, min(limits) + inset, max(limits) - inset),
  y = runif(ncircles, min(limits) + inset, max(limits) - inset),
  r = rbeta(ncircles, 1, 10) * rmax)

## ------------------------------------------------------------------------
library(packcircles)

res <- circleLayout(xyr, limits, limits, maxiter = 1000)
cat(res$niter, "iterations performed")

## ---- fig.width=7, fig.height=4------------------------------------------
library(ggplot2)
library(gridExtra)

## plot data for the `before` layout
dat.before <- circlePlotData(xyr)

## plot dta for the `after` layout returned by circleLayout
dat.after <- circlePlotData(res$layout)

doPlot <- function(dat, title)
  ggplot(dat) + 
  geom_polygon(aes(x, y, group=id), colour="brown", fill="burlywood", alpha=0.3) +
  coord_equal(xlim=limits, ylim=limits) +
  theme_bw() +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank()) +
  labs(title=title)

grid.arrange(
  doPlot(dat.before, "before"),
  doPlot(dat.after, "after"),
  nrow=1)

## ---- fig.width=7, fig.height=4------------------------------------------
largest.id <- which(xyr$r == max(xyr$r))

# add a column to the previously generated plot data for the 'before' circles
dat.before$state <- ifelse(dat.before$id == largest.id, "static", "free")

# tweak the plot function to colour circles based on the state column
doPlot <- function(dat, title)
  ggplot(dat) + 
  geom_polygon(aes(x, y, group=id, fill=state), colour="brown1") +
  scale_fill_manual(values=c("NA", "brown4")) +
  coord_equal(xlim=limits, ylim=limits) +
  theme_bw() +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        legend.position="none") +
  labs(title=title)

g.before <- doPlot(dat.before, "before")

# now re-run the layout algorithm with a weights vector to fix the position
# of the largest circle
wts <- rep(1.0, nrow(xyr))
wts[ largest.id ] <- 0.0

res <- circleLayout(xyr, limits, limits, maxiter = 1000, weights=wts)

# finally generate a new plot for the 'after' circles
dat.after <- circlePlotData(res$layout)
dat.after$state <- ifelse(dat.after$id == largest.id, "static", "free")

g.after <- doPlot(dat.after, "after")

grid.arrange(g.before, g.after, nrow=1)


