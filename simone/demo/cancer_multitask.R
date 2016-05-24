rm(list=ls())
library(simone)
data(cancer)
attach(cancer)

res.coop <- simone(expr, tasks=status)
plot(res.coop, ask=FALSE)

glist <- getNetwork(res.coop, "BIC")
plot(glist[[1]],glist[[2]])
glist <- getNetwork(res.coop, 65)
plot(glist[[1]],glist[[2]])

detach(cancer)
