## ------------------------------------------------------------------------
require(Rdistance)

## ------------------------------------------------------------------------
data(sparrow.detections)
head(sparrow.detections)

## ------------------------------------------------------------------------
sparrow.detections$dist <- perp.dists(obs.dist=sparrow.detections$sightdist,
                                      obs.angle=sparrow.detections$sightangle)
sparrow.detections <- sparrow.detections[, -which(names(sparrow.detections)
                                                  %in% c("sightdist", "sightangle"))]                                                                  
head(sparrow.detections)

## ---- fig.width=6, fig.height=4------------------------------------------
hist(sparrow.detections$dist, col="grey", main="", xlab="distance (m)")
rug(sparrow.detections$dist)
summary(sparrow.detections$dist)

## ---- fig.width=6, fig.height=4------------------------------------------
dfunc <- F.dfunc.estim(sparrow.detections, likelihood="halfnorm", w.hi=150)
plot(dfunc)
dfunc

## ------------------------------------------------------------------------
data(sparrow.transects)
head(sparrow.transects)

## ---- fig.width=6, fig.height=4------------------------------------------
fit <- F.abund.estim(dfunc, detection.data=sparrow.detections,
                     transect.data=sparrow.transects,
                     area=10000, R=100, ci=0.95, plot.bs=TRUE)
fit

## ------------------------------------------------------------------------
fit$n.hat
fit$ci

## ---- fig.width=6, fig.height=4------------------------------------------
auto <- F.automated.CDA(detection.data=sparrow.detections,
                        transect.data=sparrow.transects,
                        w.hi=150, plot=FALSE, area=10000, R=100, ci=0.95, plot.bs=TRUE)

