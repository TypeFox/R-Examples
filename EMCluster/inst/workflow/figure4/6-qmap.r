library(EMCluster)
# source("./R/q.map.new.r")
load("./data/simu.test.rda")

### quant.av
K.min <- min(ret.test[, 1]) 
K.max <- max(ret.test[, 2]) 
pv.un <- ret.test[, 14]
pv.un <- pv.un[ret.test[, 1] <= K.max & ret.test[, 2] <= K.max]

postscript("./plot/qmap.figure4.ps", width = 5, height = 5, horizontal = FALSE)
plotq(pv.un, dim = K.max - 1, start = K.min)
dev.off()
