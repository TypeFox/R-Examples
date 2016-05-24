## ---- echo = FALSE, message = FALSE--------------------------------------
##knitr::opts_chunk$set(collapse = TRUE, comment = "")
library(palr)
library(raster)
load(system.file("extdata", "oisst.rda", package = "palr"))

## ------------------------------------------------------------------------
library(raster)
library(palr)
load(system.file("extdata", "oisst.rda", package = "palr"))
plot(oisst)

## ------------------------------------------------------------------------
sstcols <- sstPal(palette = TRUE)
plot(oisst, col = sstcols$col, zlim = range(sstcols$breaks))

## ------------------------------------------------------------------------
smsst <- focal(aggregate(oisst, fact = 5, fun = mean), matrix(1, 3), fun = median, na.rm = FALSE)
psst <- rasterToPolygons(smsst, na.rm = FALSE)
layout(matrix(1:2, ncol = 2))
plot(oisst, col = sstcols$col, zlim = range(sstcols$breaks))
plot(oisst, col = "transparent", legend = FALSE)## plot twice to get exact alignment
plot(psst, col = sstPal(psst$layer), add = TRUE)
temps <- rev(c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20))
op <- par(xpd = NA); legend("right", legend = temps, fill = sstPal(temps), inset = -0.2, cex = 0.65)
par(op)

## ------------------------------------------------------------------------
metadata(oisst)$filename
writeLines(metadata(oisst)$ncdump)

