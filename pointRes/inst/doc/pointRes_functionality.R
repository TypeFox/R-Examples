## ---- echo = FALSE, message = FALSE--------------------------------------
library(pointRes)

## ------------------------------------------------------------------------
library(pointRes)
data(s033) # the result of s033 <- read.rwl('s033.rwl') - a function of the dplR package

## ------------------------------------------------------------------------
pyc <- pointer.norm(s033, window = 5, method.thresh = "Cropper", C.thresh = 0.75, series.thresh = 75)
pyn <- pointer.norm(s033, window = 5, method.thresh = "Neuwirth", N.thresh1 = 1, N.thresh2 = 1.28, N.thresh3 = 1.645, series.thresh = 75)

## ------------------------------------------------------------------------
rgc <- pointer.rgc(s033, nb.yrs = 4, rgc.thresh.pos = 60, rgc.thresh.neg = 40, series.thresh = 75)

## ------------------------------------------------------------------------
head(pyc$out) # first five lines of the summarizing 'out' component of list pyc

## ---- fig.width = 7, fig.height = 3.5, fig.retina = 3--------------------
rgc.plot(rgc, start.yr = 1950, end.yr = NULL, sd.disp = FALSE, x.tick.major = 10, x.tick.minor = 5)

## ---- fig.width = 7, fig.height = 5.3, fig.retina = 3--------------------
event.plot(pyn, start.yr = 1950, end.yr = NULL, x.tick.major = 10, x.tick.minor = 5)

## ---- fig.width = 7, fig.height = 2, fig.retina = 3----------------------
pointer.plot(list(pyn,pyn), sign = "neg", start.yr = 1950, end.yr = NULL)

## ------------------------------------------------------------------------
res <- res.comp(s033, nb.yrs = 4, res.thresh.neg = 40, series.thresh = 75)

## ------------------------------------------------------------------------
head(res$out.select)

## ---- fig.width = 6, fig.height = 5.3, fig.retina = 3--------------------
res.plot(res, select.yr = NULL, multi.panel = TRUE)

## ------------------------------------------------------------------------
citation()
citation("pointRes")


