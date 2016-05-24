library(EMCluster)
# source("./R/q.map.new.r")
load("./data/galaxy.test.rda")

### quant.av
K.min <- min(ret.test[, 1]) 
K.max <- max(ret.test[, 2]) 
pv.un <- ret.test[, 12]
pv.un <- pv.un[ret.test[, 1] <= K.max & ret.test[, 2] <= K.max]

postscript("./plot/qmap.galaxy.ps", width = 5, height = 5, horizontal = FALSE)
plotq(pv.un, dim = K.max - 1, start = K.min)
dev.off()
