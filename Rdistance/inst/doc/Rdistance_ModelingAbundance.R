## ------------------------------------------------------------------------
require(Rdistance)
data(sparrow.detections)
data(sparrow.transects)
auto <- F.automated.CDA(detection.data=sparrow.detections,
                        transect.data=sparrow.transects,
                        by.id=TRUE, w.hi=150, area=10000, R=10, ci=0.95,
                        plot=FALSE, plot.bs=FALSE)


## ------------------------------------------------------------------------
head(auto$nhat.df)

## ------------------------------------------------------------------------
mydata <- merge(sparrow.transects, auto$nhat.df, by="siteID")
head(mydata)

## ------------------------------------------------------------------------
mod <- lm(nhat ~ sagemean, data=mydata)
summary(mod)

## ---- fig.width=4.5, fig.height=5, echo=FALSE----------------------------
plot(mydata$nhat ~ mydata$sagemean,
     xlab="Sagebrush Cover (%)", ylab="Sparrow Density (birds per ha +/- 95% CI)")

pred.sage <- seq(min(mydata$sagemean), max(mydata$sagemean), by=0.5)
pred.nhat <- data.frame(predict(mod, newdata=data.frame(sagemean=pred.sage),
                                interval="confidence"))

lines(pred.sage, pred.nhat$fit, col="red", lwd=3)
lines(pred.sage, pred.nhat$lwr, col="red", lwd=2, lty="dashed")
lines(pred.sage, pred.nhat$upr, col="red", lwd=2, lty="dashed")


