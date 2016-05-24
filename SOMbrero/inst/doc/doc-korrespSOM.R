## ----loading, results='hide', echo=FALSE, warning=FALSE, message=FALSE----
library(SOMbrero)

## ----loadData------------------------------------------------------------
data(presidentielles2002)
apply(presidentielles2002, 2, sum)

## ----presiTrain, cache=TRUE----------------------------------------------
set.seed(4031719)
korresp.som <- trainSOM(x.data=presidentielles2002, dimension=c(8,8),
                        type="korresp", scaling="chi2", nb.save=10)
korresp.som

## ----energyPresi---------------------------------------------------------
plot(korresp.som, what="energy")

## ----presiClusters-------------------------------------------------------
korresp.som$clustering

## ----presiHitmap---------------------------------------------------------
plot(korresp.som, what="obs", type="hitmap")

## ----presiGraphObs, warning=FALSE, fig.width=12, fig.height=12-----------
plot(korresp.som, what="obs", type="names", scale=c(0.9,0.5))

## ----presiProtoL---------------------------------------------------------
# plot the line prototypes (106 French departements)
plot(korresp.som, what="prototypes", type="lines", view="r", print.title=TRUE)
# plot the column prototypes (16 candidates)
plot(korresp.som, what="prototypes", type="lines", view="c", print.title=TRUE)

## ----presiProtoC3d, fig.width=12, fig.height=6---------------------------
par(mfrow=c(1,2))
plot(korresp.som, what="prototypes", type="color", variable="LE_PEN")
plot(korresp.som, what="prototypes", type="3d", variable="martinique")

## ----presiProtoNumber, fig.width=12, fig.height=6------------------------
par(mfrow=c(1,2))
plot(korresp.som, what="prototypes", type="color", variable=5, view="c")
plot(korresp.som, what="prototypes", type="3d", variable=5, view="c")

## ----presiGraphProto2----------------------------------------------------
plot(korresp.som, what="prototypes", type="poly.dist", print.title=TRUE)
plot(korresp.som, what="prototypes", type="umatrix", print.title=TRUE)
plot(korresp.som, what="prototypes", type="smooth.dist", print.title=TRUE)
plot(korresp.som, what="prototypes", type="mds")
plot(korresp.som, what="prototypes", type="grid.dist")

## ----presiQuality--------------------------------------------------------
quality(korresp.som)

## ----presiSC1------------------------------------------------------------
plot(superClass(korresp.som))

## ----presiSC2------------------------------------------------------------
my.sc <- superClass(korresp.som, k=3)
summary(my.sc)
plot(my.sc, plot.var=FALSE)

## ----presiSC3------------------------------------------------------------
plot(my.sc, type="grid", plot.legend=TRUE)
plot(my.sc, type="dendro3d")

## ----presiSC4------------------------------------------------------------
plot(my.sc, type="hitmap", plot.legend=TRUE)
plot(my.sc, type="lines", print.title=TRUE)
plot(my.sc, type="lines", print.title=TRUE, view="c")
plot(my.sc, type="mds", plot.legend=TRUE)

## ----presiSC5------------------------------------------------------------
plot(my.sc, type="color", view="r", variable="correze")
plot(my.sc, type="color", view="c", variable="JOSPIN")
plot(my.sc, type="poly.dist")

## ----DOM, echo=FALSE-----------------------------------------------------
names(korresp.som$clustering)[which(korresp.som$clustering %in%
                                      which(my.sc$cluster==1))]

