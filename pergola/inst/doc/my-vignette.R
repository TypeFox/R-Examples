## ------------------------------------------------------------------------
library(pergola)

## ------------------------------------------------------------------------
set.seed(31415)

## ------------------------------------------------------------------------
data("simTetra")

## ------------------------------------------------------------------------
simTetra[1:5, 1:12]

## ------------------------------------------------------------------------
simTetraGen <- shuffleInput(simTetra, ploidy = 4)
simTetraGen[1:5, 1:12]

## ------------------------------------------------------------------------
simTetraGen <- bases2genotypes(simTetraGen, ploidy = 4)
simTetraGen[1:5, 1:6]

## ------------------------------------------------------------------------
rf <- calcRec(simTetraGen, ploidy = 4)

## ---- fig.show = 'hold', fig.width = 10, fig.height = 10-----------------
image(rf, xaxt = 'n', yaxt = 'n')
axis(side = 1, at = seq(0, 1, length.out = nrow(rf)),
     labels = rownames(rf), las = 2, cex.axis = 0.8)
axis(side = 2, at = seq(0, 1, length.out = nrow(rf)),
     labels = rownames(rf), las = 2, cex.axis = 0.8) 

## ---- fig.show = 'hold', fig.width = 15, fig.height = 5------------------
plotRf(rf)

## ---- fig.show = 'hold', fig.width = 10, fig.height = 10-----------------
plotRf(rf, plottype = "image")

## ------------------------------------------------------------------------
split <- splitChr(rf, nchr = 7)
table(split$split)
head(split)

## ---- results = 'hide'---------------------------------------------------
split <- sortLeafs(rf, split)
head(split)

## ---- fig.show = 'hold', fig.width = 10, fig.height = 10-----------------
image(rf[split$order, split$order], xaxt = 'n', yaxt = 'n')
axis(side = 1, at = seq(0, 1, length.out = nrow(rf)),
     labels = rownames(rf)[split$order], las = 2, cex.axis = 0.8)
axis(side = 2, at = seq(0, 1, length.out = nrow(rf)),
     labels = rownames(rf)[split$order], las = 2, cex.axis = 0.8) 

## ------------------------------------------------------------------------
set.seed(3)
ambRF <- cbind(c(0, 2, 4, 6, 8, 12),
               c(2, 0, 4, 4, 7, 10),
               c(4, 4, 0, 2, 4, 7),
               c(6, 4, 2, 0, 4, 5),
               c(8, 7, 4, 4, 0, 3),
               c(12, 10, 7, 5, 3, 0)) / 100
ambsplit <- data.frame(names = LETTERS[1:6],
                    split = rep(1, 6),
                    order = 1:6)
amb1 <- sortLeafs(ambRF, ambsplit)
amb1$order
amb2 <- c(1,2,4,3,5,6)
amb3 <- c(2,1,3,4,5,6)
amb4 <- 6:1 #reverse of amb1

calcSarf(ambRF, amb1$order, n = 1)
calcSarf(ambRF, amb2, n = 1)
calcSarf(ambRF, amb3, n = 1)
calcSarf(ambRF, amb4, n = 1)

calcSarf(ambRF, amb1$order, n = 2)
calcSarf(ambRF, amb2, n = 2)
calcSarf(ambRF, amb3, n = 2)
calcSarf(ambRF, amb4, n = 2)

## ------------------------------------------------------------------------
maps <- pullMap(rf, split)

## ---- fig.show='hold', fig.width = 8, fig.height = 5---------------------
plotChr(maps[[1]], cex = 0.6)

## ---- fig.show='hold', fig.width = 8, fig.height = 5---------------------
maps2 <- pullMap(rf, split, fun = "kosambi")
plotChr(maps[[1]], maps2[[1]], cex = 0.6)

## ---- results="hide",prompt=FALSE,message=FALSE--------------------------
library(dendextend)
library(gclus)

## ------------------------------------------------------------------------
maps3 <- maps2
maps3[1] <- maps2[2]
maps3[2] <- maps2[1]
maps3[3] <- maps2[7]
maps3[[4]] <- rev(max(maps3[[4]]) - maps3[[4]])
maps3[6] <- maps2[3]
maps3[7] <- maps2[6]
maps3[[4]]
maps2[[4]]

## ---- fig.show = 'hold', fig.width = 15, fig.height = 5------------------
dend1 <- map2dend(maps)
plot(dend1, cex = 0.6)

## ------------------------------------------------------------------------
dend2 <- map2dend(maps3)

## ---- fig.show = 'hold', fig.width = 10, fig.height = 15-----------------
maketangle(dend1, dend2, cutheight = 500, k = 7)

## ---- fig.show = 'hold', fig.width = 10, fig.height = 15-----------------
maps <- switchChrs(map = maps, comp = maps3)
maps <- swapChrs(map = maps, comp = maps3)
dend3 <- map2dend(maps)
dend4 <- map2dend(maps3)
maketangle(dend3, dend4, cutheight = 500, k = 7)

## ------------------------------------------------------------------------
sessionInfo()

