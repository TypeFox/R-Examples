## ----loading, results='hide', echo=FALSE, warning=FALSE, message=FALSE----
library(SOMbrero)

## ----lesmisDescr, fig.width=12, fig.height=12----------------------------
data(lesmis)
lesmis
plot(lesmis, vertex.size=0)

## ----lesmisTrain, cache=TRUE---------------------------------------------
set.seed(622)
mis.som <- trainSOM(x.data=dissim.lesmis, type="relational", nb.save=10,
                   init.proto="random", radius.type="letremy")
plot(mis.som, what="energy")

## ----lesmisClustering----------------------------------------------------
mis.som$clustering
table(mis.som$clustering)
plot(mis.som)

## ----lesmisPseudoNamesPlot, fig.height=12, fig.width=12, warning=FALSE----
plot(mis.som, what="obs", type="names")

## ----lesmisProjGraph-----------------------------------------------------
plot(mis.som, what="add", type="graph", var=lesmis)

## ----lesmisProtoProfiles-------------------------------------------------
plot(mis.som, what="prototypes", type="lines")
plot(mis.som, what="prototypes", type="radar")

## ----lesmisProtoDist-----------------------------------------------------
plot(mis.som, what="prototypes", type="poly.dist", print.title=TRUE)
plot(mis.som, what="prototypes", type="smooth.dist")
plot(mis.som, what="prototypes", type="umatrix", print.title=TRUE)
plot(mis.som, what="prototypes", type="mds")
plot(mis.som, what="prototypes", type="grid.dist")

## ----lesmisColorOverview, fig.height=12, fig.width=12--------------------
plot(lesmis, vertex.label.color=rainbow(25)[mis.som$clustering], vertex.size=0)
legend(x="left", legend=1:25, col=rainbow(25), pch=19)

## ----lesmisSCOverview----------------------------------------------------
plot(superClass(mis.som))

## ----lesmisSC------------------------------------------------------------
sc.mis <- superClass(mis.som, k=5)
summary(sc.mis)
table(sc.mis$cluster)
plot(sc.mis)
plot(sc.mis, type="grid", plot.legend=TRUE)
plot(sc.mis, type="lines", print.title=TRUE)
plot(sc.mis, type="mds", plot.legend=TRUE)
plot(sc.mis, type="dendro3d")

## ----lesmisSCColorOverview, fig.width=12, fig.height=12------------------
library(RColorBrewer)
plot(lesmis, vertex.size=0, vertex.label.color=
       brewer.pal(6, "Set2")[sc.mis$cluster[mis.som$clustering]])
legend(x="left", legend=paste("SC",1:5), col=brewer.pal(5, "Set2"), pch=19)

## ----lesmisSCProjGraph---------------------------------------------------
projectIGraph(sc.mis, lesmis)
par(mar=rep(0,4))
plot(sc.mis, type="projgraph", variable=lesmis, s.radius=2)

