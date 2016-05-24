rm(list=ls())
library(simone)

## the data set
data(cancer)
attach(cancer)

## no clustering by default
res.no <- simone(expr)
plot(res.no, ask=FALSE)
g.no <- getNetwork(res.no, 30)
plot(g.no)

## try with clustering now
control <- setOptions(clusters.crit=30)
res.cl  <- simone(expr, clustering=TRUE, control=control)
g.cl    <- getNetwork(res.cl, 30)
plot(g.cl, last.coord=TRUE)
plot(g.cl, type = "circles")

## Let us compare the two network
g.cl$name <- "clustering prior"
g.no$name <- "no clustering prior"
plot(g.no,g.cl)

detach(cancer)
