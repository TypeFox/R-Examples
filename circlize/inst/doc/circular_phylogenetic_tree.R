## ----echo = FALSE--------------------------------------------------------
library(knitr)
opts_chunk$set(fig.pos = "")

library(circlize)
circos.initialize = function(...) {
    circos.par(unit.circle.segments = 300)
    circlize::circos.initialize(...)
}

## ----eval = FALSE--------------------------------------------------------
#  library(ape)
#  data(bird.orders)
#  hc = as.hclust(bird.orders)

## ----echo = FALSE--------------------------------------------------------
load(paste0(system.file(package = "circlize"), "/extdata/bird.orders.RData"))

## ------------------------------------------------------------------------
labels = hc$labels  # name of birds
ct = cutree(hc, 6)  # cut tree into 6 pieces
n = length(labels)  # number of bird species
dend = as.dendrogram(hc)

## ----phylogenetic_tree_part1, eval = FALSE-------------------------------
#  library(circlize)
#  circos.par(cell.padding = c(0, 0, 0, 0))
#  circos.initialize(factors = "a", xlim = c(0, n)) # only one sector
#  max_height = attr(dend, "height")  # maximum height of the trees
#  circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = 0.3,
#      panel.fun = function(x, y) {
#          for(i in seq_len(n)) {
#              circos.text(i-0.5, 0, labels[i], adj = c(0, 0.5),
#                  facing = "clockwise", niceFacing = TRUE,
#                  col = ct[labels[i]], cex = 0.7)
#          }
#  })

## ----phylogenetic_tree_part2, eval = FALSE, message = FALSE, text = FALSE----
#  suppressPackageStartupMessages(library(dendextend))
#  dend = color_branches(dend, k = 6, col = 1:6)
#  circos.trackPlotRegion(ylim = c(0, max_height), bg.border = NA,
#      track.height = 0.4, panel.fun = function(x, y) {
#          circos.dendrogram(dend, max_height = max_height)
#  })
#  circos.clear()

## ----eval = FALSE--------------------------------------------------------
#  circos.dendrogram(dend, max_height = max_height, facing = "inside")

## ----phylogenetic_tree_part3, eval = FALSE, message = FALSE--------------
#  circos.par(cell.padding = c(0, 0, 0, 0))
#  circos.initialize(factors = "a", xlim = c(0, n)) # only one sector
#  circos.trackPlotRegion(ylim = c(0, max_height), bg.border = NA,
#      track.height = 0.4, panel.fun = function(x, y) {
#          circos.dendrogram(dend, max_height = max_height, facing = "inside")
#  })
#  circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = 0.3,
#      panel.fun = function(x, y) {
#          for(i in seq_len(n)) {
#              circos.text(i-0.5, 1, labels[i], adj = c(1, 0.5),
#                  facing = "clockwise", niceFacing = TRUE,
#                  col = ct[labels[i]], cex = 0.7)
#          }
#  })
#  circos.clear()

## ----phylogenetic_tree, echo = FALSE, fig.align = "center", out.width = "\\textwidth", out.height = "\\textwidth", fig.width = 10, fig.height = 10, fig.cap = "A simple phylogenetic tree. A: circular style; B) dendrogram is facing inside of the circle; C: normal style."----
par(mfrow = c(2, 2))
library(circlize)
circos.par(cell.padding = c(0, 0, 0, 0))
circos.initialize(factors = "a", xlim = c(0, n)) # only one sector
max_height = attr(dend, "height")  # maximum height of the trees
circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = 0.3, 
    panel.fun = function(x, y) {
        for(i in seq_len(n)) {
            circos.text(i-0.5, 0, labels[i], adj = c(0, 0.5), 
                facing = "clockwise", niceFacing = TRUE,
                col = ct[labels[i]], cex = 0.7)
        }
})
suppressPackageStartupMessages(library(dendextend))
dend = color_branches(dend, k = 6, col = 1:6)
circos.trackPlotRegion(ylim = c(0, max_height), bg.border = NA, 
    track.height = 0.4, panel.fun = function(x, y) {
        circos.dendrogram(dend, max_height = max_height)
})
circos.clear()
text(-0.9, 0.9, "A", cex = 1.5)
circos.par(cell.padding = c(0, 0, 0, 0))
circos.initialize(factors = "a", xlim = c(0, n)) # only one sector
circos.trackPlotRegion(ylim = c(0, max_height), bg.border = NA, 
    track.height = 0.4, panel.fun = function(x, y) {
        circos.dendrogram(dend, max_height = max_height, facing = "inside")
})
circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = 0.3, 
    panel.fun = function(x, y) {
        for(i in seq_len(n)) {
            circos.text(i-0.5, 1, labels[i], adj = c(1, 0.5), 
                facing = "clockwise", niceFacing = TRUE,
                col = ct[labels[i]], cex = 0.7)
        }
})
circos.clear()
text(-0.9, 0.9, "B", cex = 1.5)

par(mar = c(8, 4, 4, 1))
plot(dend, main = "C")

