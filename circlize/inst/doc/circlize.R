## ----echo = FALSE--------------------------------------------------------
library(knitr)
opts_chunk$set(fig.pos = "", fig.align = "center")

library(circlize)
circos.initialize = function(...) {
    circos.par(unit.circle.segments = 200)
    circlize::circos.initialize(...)
}

## ----circlize_correspondence, echo = FALSE, fig.cap = "Correspondence between graphic functions in \\textbf{circlize} and in traditional R graphic engine."----
source("src/intro-00-correspondence.R")

## ------------------------------------------------------------------------
set.seed(999)
n = 1000
a = data.frame(factor = sample(letters[1:8], n, replace = TRUE),
    x = rnorm(n), y = runif(n))

## ----circlize_glance_0, eval = FALSE-------------------------------------
#  library(circlize)
#  par(mar = c(1, 1, 1, 1), lwd = 0.1, cex = 0.7)
#  circos.par("track.height" = 0.1)
#  circos.initialize(factors = a$factor, x = a$x)

## ----circlize_glance_1, eval=FALSE---------------------------------------
#  circos.trackPlotRegion(factors = a$factor, y = a$y,
#      panel.fun = function(x, y) {
#          circos.axis()
#  })
#  col = rep(c("#FF0000", "#00FF00"), 4)
#  circos.trackPoints(a$factor, a$x, a$y, col = col, pch = 16, cex = 0.5)
#  circos.text(-1, 0.5, "left", sector.index = "a", track.index = 1)
#  circos.text(1, 0.5, "right", sector.index = "a")

## ----circlize_glance_2, eval=FALSE---------------------------------------
#  bgcol = rep(c("#EFEFEF", "#CCCCCC"), 4)
#  circos.trackHist(a$factor, a$x, bg.col = bgcol, col = NA)

## ----circlize_glance_3, eval=FALSE---------------------------------------
#  circos.trackPlotRegion(factors = a$factor, x = a$x, y = a$y,
#      panel.fun = function(x, y) {
#          grey = c("#FFFFFF", "#CCCCCC", "#999999")
#          sector.index = get.cell.meta.data("sector.index")
#          xlim = get.cell.meta.data("xlim")
#          ylim = get.cell.meta.data("ylim")
#          circos.text(mean(xlim), mean(ylim), sector.index)
#          circos.points(x[1:10], y[1:10], col = "red", pch = 16, cex = 0.6)
#          circos.points(x[11:20], y[11:20], col = "blue", cex = 0.6)
#  })

## ----circlize_glance_3_update, eval=FALSE--------------------------------
#  circos.updatePlotRegion(sector.index = "d", track.index = 2)
#  circos.points(x = -2:2, y = rep(0, 5))
#  xlim = get.cell.meta.data("xlim")
#  ylim = get.cell.meta.data("ylim")
#  circos.text(mean(xlim), mean(ylim), "updated")

## ----circlize_glance_4, eval=FALSE---------------------------------------
#  circos.trackPlotRegion(factors = a$factor, y = a$y)
#  circos.trackLines(a$factor[1:100], a$x[1:100], a$y[1:100], type = "h")

## ----circlize_glance_5, eval=FALSE---------------------------------------
#  circos.link("a", 0, "b", 0, h = 0.4)
#  circos.link("c", c(-0.5, 0.5), "d", c(-0.5,0.5), col = "red",
#      border = "blue", h = 0.2)
#  circos.link("e", 0, "g", c(-1,1), col = "green", border = "black", lwd = 2, lty = 2)

## ----circlize_glance, echo = FALSE, out.width = "0.6\\textheight", out.height = "0.9\\textheight", fig.width = 7, fig.height = 10.5, fig.cap = "A step-by-step example by \\textbf{circlize}."----
par(mfrow = c(3, 2), mar = c(1, 1, 1, 1))

library(circlize)
par(mar = c(1, 1, 1, 1), lwd = 0.1, cex = 0.7)
circos.par("track.height" = 0.1)
circos.initialize(factors = a$factor, x = a$x)
circos.trackPlotRegion(factors = a$factor, y = a$y,
    panel.fun = function(x, y) {
        circos.axis()
})
col = rep(c("#FF0000", "#00FF00"), 4)
circos.trackPoints(a$factor, a$x, a$y, col = col, pch = 16, cex = 0.5)
circos.text(-1, 0.5, "left", sector.index = "a", track.index = 1)
circos.text(1, 0.5, "right", sector.index = "a")
circos.clear()
text(-0.9, 0.9, "A", cex = 1.5)

library(circlize)
par(mar = c(1, 1, 1, 1), lwd = 0.1, cex = 0.7)
circos.par("track.height" = 0.1)
circos.initialize(factors = a$factor, x = a$x)
circos.trackPlotRegion(factors = a$factor, y = a$y,
    panel.fun = function(x, y) {
        circos.axis()
})
col = rep(c("#FF0000", "#00FF00"), 4)
circos.trackPoints(a$factor, a$x, a$y, col = col, pch = 16, cex = 0.5)
circos.text(-1, 0.5, "left", sector.index = "a", track.index = 1)
circos.text(1, 0.5, "right", sector.index = "a")
bgcol = rep(c("#EFEFEF", "#CCCCCC"), 4)
circos.trackHist(a$factor, a$x, bg.col = bgcol, col = NA)
circos.clear()
text(-0.9, 0.9, "B", cex = 1.5)

library(circlize)
par(mar = c(1, 1, 1, 1), lwd = 0.1, cex = 0.7)
circos.par("track.height" = 0.1)
circos.initialize(factors = a$factor, x = a$x)
circos.trackPlotRegion(factors = a$factor, y = a$y,
    panel.fun = function(x, y) {
        circos.axis()
})
col = rep(c("#FF0000", "#00FF00"), 4)
circos.trackPoints(a$factor, a$x, a$y, col = col, pch = 16, cex = 0.5)
circos.text(-1, 0.5, "left", sector.index = "a", track.index = 1)
circos.text(1, 0.5, "right", sector.index = "a")
bgcol = rep(c("#EFEFEF", "#CCCCCC"), 4)
circos.trackHist(a$factor, a$x, bg.col = bgcol, col = NA)
circos.trackPlotRegion(factors = a$factor, x = a$x, y = a$y,
    panel.fun = function(x, y) {
        grey = c("#FFFFFF", "#CCCCCC", "#999999")
        sector.index = get.cell.meta.data("sector.index")
        xlim = get.cell.meta.data("xlim")
        ylim = get.cell.meta.data("ylim")
        circos.text(mean(xlim), mean(ylim), sector.index)
        circos.points(x[1:10], y[1:10], col = "red", pch = 16, cex = 0.6)
        circos.points(x[11:20], y[11:20], col = "blue", cex = 0.6)
})
circos.clear()
text(-0.9, 0.9, "C", cex = 1.5)

library(circlize)
par(mar = c(1, 1, 1, 1), lwd = 0.1, cex = 0.7)
circos.par("track.height" = 0.1)
circos.initialize(factors = a$factor, x = a$x)
circos.trackPlotRegion(factors = a$factor, y = a$y,
    panel.fun = function(x, y) {
        circos.axis()
})
col = rep(c("#FF0000", "#00FF00"), 4)
circos.trackPoints(a$factor, a$x, a$y, col = col, pch = 16, cex = 0.5)
circos.text(-1, 0.5, "left", sector.index = "a", track.index = 1)
circos.text(1, 0.5, "right", sector.index = "a")
bgcol = rep(c("#EFEFEF", "#CCCCCC"), 4)
circos.trackHist(a$factor, a$x, bg.col = bgcol, col = NA)
circos.trackPlotRegion(factors = a$factor, x = a$x, y = a$y,
    panel.fun = function(x, y) {
        grey = c("#FFFFFF", "#CCCCCC", "#999999")
        sector.index = get.cell.meta.data("sector.index")
        xlim = get.cell.meta.data("xlim")
        ylim = get.cell.meta.data("ylim")
        circos.text(mean(xlim), mean(ylim), sector.index)
        circos.points(x[1:10], y[1:10], col = "red", pch = 16, cex = 0.6)
        circos.points(x[11:20], y[11:20], col = "blue", cex = 0.6)
})
circos.updatePlotRegion(sector.index = "d", track.index = 2)
circos.points(x = -2:2, y = rep(0, 5))
xlim = get.cell.meta.data("xlim")
ylim = get.cell.meta.data("ylim")
circos.text(mean(xlim), mean(ylim), "updated")
circos.clear()
text(-0.9, 0.9, "D", cex = 1.5)

library(circlize)
par(mar = c(1, 1, 1, 1), lwd = 0.1, cex = 0.7)
circos.par("track.height" = 0.1)
circos.initialize(factors = a$factor, x = a$x)
circos.trackPlotRegion(factors = a$factor, y = a$y,
    panel.fun = function(x, y) {
        circos.axis()
})
col = rep(c("#FF0000", "#00FF00"), 4)
circos.trackPoints(a$factor, a$x, a$y, col = col, pch = 16, cex = 0.5)
circos.text(-1, 0.5, "left", sector.index = "a", track.index = 1)
circos.text(1, 0.5, "right", sector.index = "a")
bgcol = rep(c("#EFEFEF", "#CCCCCC"), 4)
circos.trackHist(a$factor, a$x, bg.col = bgcol, col = NA)
circos.trackPlotRegion(factors = a$factor, x = a$x, y = a$y,
    panel.fun = function(x, y) {
        grey = c("#FFFFFF", "#CCCCCC", "#999999")
        sector.index = get.cell.meta.data("sector.index")
        xlim = get.cell.meta.data("xlim")
        ylim = get.cell.meta.data("ylim")
        circos.text(mean(xlim), mean(ylim), sector.index)
        circos.points(x[1:10], y[1:10], col = "red", pch = 16, cex = 0.6)
        circos.points(x[11:20], y[11:20], col = "blue", cex = 0.6)
})
circos.updatePlotRegion(sector.index = "d", track.index = 2)
circos.points(x = -2:2, y = rep(0, 5))
xlim = get.cell.meta.data("xlim")
ylim = get.cell.meta.data("ylim")
circos.text(mean(xlim), mean(ylim), "updated")
circos.trackPlotRegion(factors = a$factor, y = a$y)
circos.trackLines(a$factor[1:100], a$x[1:100], a$y[1:100], type = "h")
circos.clear()
text(-0.9, 0.9, "E", cex = 1.5)

library(circlize)
par(mar = c(1, 1, 1, 1), lwd = 0.1, cex = 0.7)
circos.par("track.height" = 0.1)
circos.initialize(factors = a$factor, x = a$x)
circos.trackPlotRegion(factors = a$factor, y = a$y,
    panel.fun = function(x, y) {
        circos.axis()
})
col = rep(c("#FF0000", "#00FF00"), 4)
circos.trackPoints(a$factor, a$x, a$y, col = col, pch = 16, cex = 0.5)
circos.text(-1, 0.5, "left", sector.index = "a", track.index = 1)
circos.text(1, 0.5, "right", sector.index = "a")
bgcol = rep(c("#EFEFEF", "#CCCCCC"), 4)
circos.trackHist(a$factor, a$x, bg.col = bgcol, col = NA)
circos.trackPlotRegion(factors = a$factor, x = a$x, y = a$y,
    panel.fun = function(x, y) {
        grey = c("#FFFFFF", "#CCCCCC", "#999999")
        sector.index = get.cell.meta.data("sector.index")
        xlim = get.cell.meta.data("xlim")
        ylim = get.cell.meta.data("ylim")
        circos.text(mean(xlim), mean(ylim), sector.index)
        circos.points(x[1:10], y[1:10], col = "red", pch = 16, cex = 0.6)
        circos.points(x[11:20], y[11:20], col = "blue", cex = 0.6)
})
circos.updatePlotRegion(sector.index = "d", track.index = 2)
circos.points(x = -2:2, y = rep(0, 5))
xlim = get.cell.meta.data("xlim")
ylim = get.cell.meta.data("ylim")
circos.text(mean(xlim), mean(ylim), "updated")
circos.trackPlotRegion(factors = a$factor, y = a$y)
circos.trackLines(a$factor[1:100], a$x[1:100], a$y[1:100], type = "h")
circos.link("a", 0, "b", 0, h = 0.4)
circos.link("c", c(-0.5, 0.5), "d", c(-0.5,0.5), col = "red",
    border = "blue", h = 0.2)
circos.link("e", 0, "g", c(-1,1), col = "green", border = "black", lwd = 2, lty = 2)
circos.clear()
text(-0.9, 0.9, "F", cex = 1.5)
par(mfrow = c(1, 1))

## ----eval = FALSE--------------------------------------------------------
#  circos.info()
#  circos.info(sector.index = "a", track.index = 2)

## ----eval=FALSE----------------------------------------------------------
#  circos.clear()

## ----circlize_coordinate_transformation, echo = FALSE, out.width = "0.3546099\\textheight", out.height = "0.8\\textheight", fig.width = 5, fig.height = 11.28, fig.cap = "Transformation between different coordinates. Top: data coordinate; Middle: polar coordinate; Bottom: canvas coordinate."----
source("src/intro-03-transformation.R")

## ----circlize_order, echo = FALSE, out.width = "\\textwidth", out.height = "0.8571429\\textwidth", fig.height = 6, fig.width = 7, fig.cap = "Order of drawing circular layout."----
source("src/intro-02-order.R")

## ----eval = FALSE--------------------------------------------------------
#  circos.initialize(factors, xlim)
#  circos.trackPlotRegion(factors, ylim)
#  for(sector.index in all.sector.index) {
#      circos.points(x1, y1, sector.index)
#      circos.lines(x2, y2, sector.index)
#  }

## ----eval = FALSE--------------------------------------------------------
#  circos.initialize(factors, xlim)
#  circos.trackPlotRegion(factors, ylim)
#  circos.trackPoints(factors, x, y)
#  circos.trackLines(factors, x, y)

## ----eval = FALSE--------------------------------------------------------
#  circos.initialize(factors, xlim)
#  circos.trackPlotRegion(factors, all_x, all_y, ylim,
#      panel.fun = function(x, y) {
#          circos.points(x, y)
#          circos.lines(x, y)
#  })

## ----circlize_coordinate, echo = FALSE, out.width = "0.8\\textwidth", fig.cap = "Sectors and tracks in circular layout. There are 10 sectors and 4 tracks. Orders of sectors are randomly permuted."----
source("src/intro-04-coordinate.R")

## ----eval = FALSE--------------------------------------------------------
#  circos.initialize(factors, x)
#  circos.initialize(factors, xlim)

## ----eval = FALSE--------------------------------------------------------
#  fa = c("d", "f", "e", "c", "g", "b", "a")
#  f1 = factor(fa)
#  circos.initialize(factors = f1, xlim = c(0, 1))
#  f2 = factor(fa, levels = fa)
#  circos.initialize(factors = f2, xlim = c(0, 1))

## ----circlize_factor, echo = FALSE, out.width = "\\textwidth", out.height = "0.5\\textwidth", fig.width = 10, fig.height = 5, fig.cap = "Different sector orders."----
source("src/intro-05-factor.R")

## ----eval = FALSE--------------------------------------------------------
#  circos.trackPlotRegion(factors, y)
#  circos.trackPlotRegion(factors, ylim)

## ----circlize_region, echo = FALSE, out.width = "0.8\\textwidth", fig.cap = "Regions in a cell"----
source("src/intro-06-region.R")

## ----circlize_direction, echo = FALSE, out.width = "\\textwidth", out.height = "0.5\\textwidth", fig.width = 10, fig.height = 5, fig.cap = "Sector directions. Sector orders are {\\tt a, b, ..., h}."----
source("src/intro-07-direction.R")

## ----eval = FALSE--------------------------------------------------------
#  circos.trackPlotRegion(data, ylim = c(0, 1), track.index = 1, ...)

## ----eval = FALSE--------------------------------------------------------
#  circos.updatePlotRegion(sector.index, track.index)
#  circos.points(x, y, sector.index, track.index)

## ----eval = FALSE--------------------------------------------------------
#  circos.points(x, y)
#  circos.points(x, y, sector.index, track.index)
#  circos.points(x, y, pch, col, cex)

## ----circlize_lines, echo = FALSE, out.width = "0.8\\textwidth", fig.cap = "Line styles"----
source("src/intro-08-lines.R")

## ----circlize_linecurve, echo = FALSE, out.width = "0.8\\textwidth", fig.cap = "Straight lines will be transformed into curves in the circle."----
source("src/intro-08-linescurve.R")

## ----eval = FALSE--------------------------------------------------------
#  circos.lines(x, y)
#  circos.lines(x, y, sector.index, track.index)
#  circos.lines(x, y, col, lwd, lty, type, straight)
#  circos.lines(x, y, col, area, baseline, border)

## ----circlize_text, echo = FALSE, out.width = "0.8\\textwidth", fig.cap = "Text facing."----
source("src/intro-09-text.R")

## ----eval = FALSE--------------------------------------------------------
#  circos.text(x, y, labels)
#  circos.text(x, y, labels, sector.index, track.index)
#  circos.text(x, y, labels, facing, adj, cex, col, font)

## ----eval = FALSE--------------------------------------------------------
#  circos.text(x, y, labels, adj = c(0, degree(5)), facing = "clockwise")

## ----circlize_text_easy, echo = FALSE, out.width = "0.6\\textheight", out.height = "0.9\\textheight", fig.width = 6, fig.height = 9, fig.cap = "Human easy text facing. When {\\tt niceFacing} is on, settings in the same row are actually identical. Red dots represent positions of the texts."----
source("src/intro-09-text-niceFacing.R")

## ----eval = FALSE--------------------------------------------------------
#  circos.rect(xleft, ybottom, xright, ytop)
#  circos.rect(xleft, ybottom, xright, ytop, sector.index, track.index)
#  circos.rect(xleft, ybottom, xright, ytop, col, border, lty, lwd)

## ----eval = FALSE--------------------------------------------------------
#  circos.polygon(x, y)
#  circos.polygon(x, y, sector.index, track.index)
#  circos.polygon(x, y, col, border, lty, lwd)

## ----circlize_errorline, echo = FALSE, out.width = "0.8\\textwidth", fig.cap = "Area of standard deviation of the smoothed line."----
source("src/intro-10-smooth.R")

## ----circlize_axis, echo = FALSE, out.width = "0.6\\textwidth", out.height = "1.2\\textwidth", fig.width = 6, fig.height = 12, fig.cap = "Axes"----
source("src/intro-11-axis.R")

## ----eval = FALSE--------------------------------------------------------
#  circos.axis(h)
#  circos.axis(h, sector.index, track.index)
#  circos.axis(h, major.at, labels, major.tick)
#  circos.axis(h, major.at, labels, major.tick, labels.font, labels.cex,
#              labels.facing, labels.away.percentage)
#  circos.axis(h, major.at, labels, major.tick, minor.ticks,
#              major.tick.percentage, lwd)

## ----eval = FALSE--------------------------------------------------------
#  circos.yaxis(side)
#  circos.yaxis(at, labels, sector.index, track.index)

## ----eval = FALSE--------------------------------------------------------
#  circos.link(sector.index1, 0, sector.index2, 0)
#  circos.link(sector.index1, c(0, 1), sector.index2, 0)
#  circos.link(sector.index1, c(0, 1), sector.index2, c(1, 2))
#  circos.link(sector.index1, c(0, 1), sector.index2, 0, col, lwd, lty, border)

## ----circlize_link, echo = FALSE, out.width = "0.6\\textheight", out.height = "0.9\\textheight", fig.width = 6, fig.height = 9, fig.cap = "A) set different positions of roots; B) set different height of two borders. C,D) set different {\\tt h} and {\\tt w}. E) if two branches overlap, the link will be degenerated as a 'hill'. F) links with directions."----
source("src/intro-12-link.R")

## ----tidy = TRUE---------------------------------------------------------
circlize:::get_most_inside_radius

## ----eval = FALSE--------------------------------------------------------
#  circos.link(sector.index1, 0, sector.index2, 0, rou)
#  circos.link(sector.index1, 0, sector.index2, 0, rou1, rou2)

## ----eval = FALSE--------------------------------------------------------
#  circos.link(sector.index1, 0, sector.index2, 0, h)
#  circos.link(sector.index1, 0, sector.index2, 0, h, h2)

## ----eval = FALSE--------------------------------------------------------
#  circos.link(sector.index1, 0, sector.index2, 0, w)
#  circos.link(sector.index1, 0, sector.index2, 0, w, w2)

## ----eval = FALSE--------------------------------------------------------
#  circos.link(sector.index1, 0, sector.index2, 0, directional = 1)
#  circos.link(sector.index1, c(0, 1), sector.index2, c(0, 1), directional = -1)

## ----eval = FALSE--------------------------------------------------------
#  factors = c("a", "a", "a", "b", "b")
#  x = 1:5
#  y = 5:1
#  circos.trackPlotRegion(factors = factors, x = x, y = y,
#      panel.fun = function(x, y) {
#          circos.points(x, y)
#  })

## ----circlize_hist, echo = FALSE, out.width = "0.8\\textwidth", fig.cap = "Histograms in circular layout."----
source("src/intro-14-hist.R")

## ----circlize_heatmap, echo = FALSE, out.width = "0.8\\textwidth", fig.cap = "Circular heatmap with dendrogram trees."----
source("src/intro-15-heatmap.R")

## ----eval = FALSE--------------------------------------------------------
#  get.cell.meta.data(name)
#  get.cell.meta.data(name, sector.index, track.index)

## ----echo = 2:8----------------------------------------------------------
pdf(NULL)
factors = c("a", "b")
circos.initialize(factors, xlim = c(0, 1))
circos.trackPlotRegion(ylim = c(0, 1))
circlize(0.5, 0.5, sector.index = "a", track.index = 1)
reverse.circlize(90, 0.9, sector.index = "a", track.index = 1)
reverse.circlize(90, 0.9, sector.index = "b", track.index = 1)
circos.clear()
invisible(dev.off())

## ----eval = FALSE--------------------------------------------------------
#  draw.sector(start.degree, end.degree, rou1)
#  draw.sector(start.degree, end.degree, rou1, rou2, center)
#  draw.sector(start.degree, end.degree, rou1, rou2, center, col, border, lwd, lty)

## ----eval = FALSE--------------------------------------------------------
#  draw.sector(start.degree, end.degree, clock.wise = FALSE)

## ----circlize_draw_sector_general, out.width = "0.8\\textwidth", fig.cap = "Examples of {\\tt draw.sector}."----
par(mar = c(1, 1, 1, 1))
plot(c(-1, 1), c(-1, 1), type = "n", axes = FALSE, ann = FALSE)
draw.sector(20, 0)
draw.sector(30, 60, rou1 = 0.8, rou2 = 0.5, clock.wise = FALSE, col = "#FF000080")
draw.sector(350, 1000, col = "#00FF0080", border = NA)
draw.sector(0, 180, rou1 = 0.25, center = c(-0.5, 0.5), border = 2, lwd = 2, lty = 2)
draw.sector(0, 360, rou1 = 0.7, rou2 = 0.6, col = "#0000FF80")

## ----circlize_highlight_1, eval = FALSE----------------------------------
#  par(mar = c(1, 1, 1, 1))
#  factors = letters[1:8]
#  circos.initialize(factors, xlim = c(0, 1))
#  for(i in 1:3) {
#      circos.trackPlotRegion(ylim = c(0, 1))
#  }
#  circos.info(plot = TRUE)

## ----circlize_highlight_2, eval = FALSE----------------------------------
#  draw.sector(get.cell.meta.data("cell.start.degree", sector.index = "a"),
#              get.cell.meta.data("cell.end.degree", sector.index = "a"),
#              rou1 = 1, col = "#FF000040")

## ----circlize_highlight_3, eval = FALSE----------------------------------
#  draw.sector(0, 360,
#      rou1 = get.cell.meta.data("cell.top.radius", track.index = 1),
#      rou2 = get.cell.meta.data("cell.bottom.radius", track.index = 1),
#      col = "#00FF0040")

## ----circlize_highlight_4, eval = FALSE----------------------------------
#  draw.sector(get.cell.meta.data("cell.start.degree", sector.index = "e"),
#              get.cell.meta.data("cell.end.degree", sector.index = "f"),
#              get.cell.meta.data("cell.top.radius", track.index = 2),
#              get.cell.meta.data("cell.bottom.radius", track.index = 3),
#              col = "#0000FF40")

## ----circlize_highlight_5, eval = FALSE----------------------------------
#  pos = circlize(c(0.2, 0.8), c(0.2, 0.8), sector.index = "h", track.index = 2)
#  draw.sector(pos[1, "theta"], pos[2, "theta"], pos[1, "rou"], pos[2, "rou"],
#      clock.wise = TRUE, col = "#00FFFF40")
#  circos.clear()

## ----circlize_highlight_sector, eval = FALSE-----------------------------
#  factors = letters[1:8]
#  circos.initialize(factors, xlim = c(0, 1))
#  for(i in 1:4) {
#      circos.trackPlotRegion(ylim = c(0, 1))
#  }
#  circos.info(plot = TRUE)
#  
#  highlight.sector(c("a", "h"), track.index = 1, text = "a and h belong to a same group",
#      facing = "bending.inside", niceFacing = TRUE, text.vjust = -3)
#  highlight.sector("c", col = "#00FF0040")
#  highlight.sector("d", col = NA, border = "red", lwd = 2)
#  highlight.sector("e", col = "#0000FF40", track.index = c(2, 3))
#  highlight.sector(c("f", "g"), col = NA, border = "green",
#      lwd = 2, track.index = c(2, 3))
#  highlight.sector(factors, col = "#FFFF0040", track.index = 4)
#  circos.clear()

## ----circlize_highlight, echo = FALSE, out.width = "0.6\\textwidth", out.height = "1.2\\textwidth", fig.width = 6, fig.height = 12, fig.cap = "Highlight sectors and tracks. A) highlight by {\\tt code draw.sector}; B) highlight by {\\tt highlight.sector}."----
par(mfrow = c(2, 1))
par(mar = c(1, 1, 1, 1))
factors = letters[1:8]
circos.initialize(factors, xlim = c(0, 1))
for(i in 1:3) {
    circos.trackPlotRegion(ylim = c(0, 1))
}
circos.info(plot = TRUE)
draw.sector(get.cell.meta.data("cell.start.degree", sector.index = "a"),
            get.cell.meta.data("cell.end.degree", sector.index = "a"),
            rou1 = 1, col = "#FF000040")
draw.sector(0, 360, 
    rou1 = get.cell.meta.data("cell.top.radius", track.index = 1),
    rou2 = get.cell.meta.data("cell.bottom.radius", track.index = 1),
    col = "#00FF0040")           
draw.sector(get.cell.meta.data("cell.start.degree", sector.index = "e"),
            get.cell.meta.data("cell.end.degree", sector.index = "f"),
            get.cell.meta.data("cell.top.radius", track.index = 2),
            get.cell.meta.data("cell.bottom.radius", track.index = 3),
            col = "#0000FF40")
pos = circlize(c(0.2, 0.8), c(0.2, 0.8), sector.index = "h", track.index = 2)
draw.sector(pos[1, "theta"], pos[2, "theta"], pos[1, "rou"], pos[2, "rou"], 
    clock.wise = TRUE, col = "#00FFFF40")
circos.clear()
text(-0.9, 0.9, "A", cex = 1.5)
factors = letters[1:8]
circos.initialize(factors, xlim = c(0, 1))
for(i in 1:4) {
    circos.trackPlotRegion(ylim = c(0, 1))
}
circos.info(plot = TRUE)

highlight.sector(c("a", "h"), track.index = 1, text = "a and h belong to a same group",
    facing = "bending.inside", niceFacing = TRUE, text.vjust = -3)
highlight.sector("c", col = "#00FF0040")
highlight.sector("d", col = NA, border = "red", lwd = 2)
highlight.sector("e", col = "#0000FF40", track.index = c(2, 3))
highlight.sector(c("f", "g"), col = NA, border = "green", 
    lwd = 2, track.index = c(2, 3))
highlight.sector(factors, col = "#FFFF0040", track.index = 4)
circos.clear()
text(-0.9, 0.9, "B", cex = 1.5)
par(mfrow = c(1, 1))

## ----echo = 2:7----------------------------------------------------------
pdf(NULL)
factors = letters[1:3]
circos.initialize(factors = factors, xlim = c(1, 2))
circos.info()
circos.trackPlotRegion(ylim = c(0, 1))
circos.info(sector.index = "a", track.index = 1)
circos.clear()
invisible(dev.off())

## ----circlize_combine, out.width = "0.8\\textwidth", fig.cap = "Combine low-level graphic functions to generate high-level graphics."----
category = paste0("category", "_", 1:10)
percent = sort(sample(40:80, 10))
color = rev(rainbow(length(percent)))

par(mar = c(1, 1, 1, 1))
circos.par("start.degree" = 90)
circos.initialize("a", xlim = c(0, 100)) # 'a` just means there is one sector
circos.trackPlotRegion(ylim = c(0.5, length(percent)+0.5), , track.height = 0.8, 
    bg.border = NA, panel.fun = function(x, y) {
        xlim = get.cell.meta.data("xlim") # in fact, it is c(0, 100)
        for(i in seq_along(percent)) {
            circos.lines(xlim, c(i, i), col = "#CCCCCC")
            circos.rect(0, i - 0.45, percent[i], i + 0.45, col = color[i], 
                border = "white")
        }

        for(i in seq_along(percent)) {
            circos.text(xlim[1], i, paste0(category[i], " - ", percent[i], "%"), 
                facing = "downward", adj = c(1.1, 0.5)) 
        }

        breaks = seq(0, 90, by = 5)
        circos.axis(h = "top", major.at = breaks, labels = paste0(breaks, "%"),
            major.tick.percentage = 0.02, labels.cex = 0.6, 
                labels.away.percentage = 0.01)
})
circos.clear()

## ------------------------------------------------------------------------
df = data.frame(factors = sample(letters[1:6], 100, replace = TRUE),
                x = rnorm(100),
                y = rnorm(100),
                stringsAsFactors = FALSE)

## ------------------------------------------------------------------------
zoom_df = df[df$factors %in% c("a", "b"), ]

## ------------------------------------------------------------------------
zoom_df$factors = paste0("zoom_", zoom_df$factors)
df2 = rbind(df, zoom_df)

## ------------------------------------------------------------------------
xrange = tapply(df2$x, df2$factors, function(x) max(x) - min(x))
normal_sector_index = unique(df$factors)
zoomed_sector_index = unique(zoom_df$factors)
sector.width = c(xrange[normal_sector_index] / sum(xrange[normal_sector_index]), 
                 xrange[zoomed_sector_index] / sum(xrange[zoomed_sector_index]))

## ----circlize_zoom_1, eval = FALSE---------------------------------------
#  par(mar = c(1, 1, 1, 1))
#  circos.par(start.degree = 90)
#  circos.initialize(df2$factors, x = df2$x, sector.width = sector.width)
#  circos.trackPlotRegion(df2$factors, x = df2$x, y = df2$y, panel.fun = function(x, y) {
#      circos.points(x, y, col = "red", pch = 16, cex = 0.5)
#      xlim = get.cell.meta.data("xlim")
#      ylim = get.cell.meta.data("ylim")
#      sector.index = get.cell.meta.data("sector.index")
#      circos.text(mean(xlim), mean(ylim), sector.index, niceFacing = TRUE)
#  })

## ----circlize_zoom_2, eval = FALSE---------------------------------------
#  circos.link("a", get.cell.meta.data("cell.xlim", sector.index = "a"),
#      "zoom_a", get.cell.meta.data("cell.xlim", sector.index = "zoom_a"),
#      border = NA, col = "#00000020")
#  circos.clear()

## ----circlize_zoom, echo = FALSE, out.width = "0.8\\textwidth", fig.cap = "Zoom sectors."----
par(mar = c(1, 1, 1, 1))
circos.par(start.degree = 90)
circos.initialize(df2$factors, x = df2$x, sector.width = sector.width)
circos.trackPlotRegion(df2$factors, x = df2$x, y = df2$y, panel.fun = function(x, y) {
    circos.points(x, y, col = "red", pch = 16, cex = 0.5)
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.index = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), mean(ylim), sector.index, niceFacing = TRUE)
})
circos.link("a", get.cell.meta.data("cell.xlim", sector.index = "a"),
    "zoom_a", get.cell.meta.data("cell.xlim", sector.index = "zoom_a"),
    border = NA, col = "#00000020")
circos.clear()

## ----eval = FALSE--------------------------------------------------------
#  par(mar = c(1, 1, 1, 1))
#  circos.par("canvas.xlim" = c(0, 1), "canvas.ylim" = c(0, 1),
#      "clock.wise" = FALSE, "gap.degree" = 0)
#  factors = letters[1:4]
#  circos.initialize(factors = factors, xlim = c(0, 1))
#  circos.trackPlotRegion(factors = factors, ylim = c(0, 1), bg.border = NA)
#  circos.updatePlotRegion(sector.index = "a", bg.border = "black")
#  x1 = runif(100)
#  y1 = runif(100)
#  circos.points(x1, y1, pch = 16, cex = 0.5)
#  circos.trackPlotRegion(factors = factors, ylim = c(0, 1), bg.border = NA)
#  circos.updatePlotRegion(sector.index = "a", bg.border = "black")
#  circos.lines(1:100/100, y1, pch = 16, cex = 0.5)
#  circos.clear()

## ----circlize_part, echo = FALSE, out.width = "0.6\\textwidth", out.height = "1.2\\textwidth", fig.width = 6, fig.height = 12, fig.cap = "One quarter of the circle."----
source("src/intro-17-part.R")

## ----eval=FALSE----------------------------------------------------------
#  par(mar = c(1, 1, 1, 1))
#  factors = letters[1:4]
#  circos.initialize(factors = factors, xlim = c(0, 1))
#  circos.trackPlotRegion(factors = factors, ylim = c(0, 1), bg.col = NA, bg.border = NA)
#  circos.updatePlotRegion(sector.index = "a", bg.border = "black")
#  x1 = runif(100)
#  y1 = runif(100)
#  circos.points(x1, y1, pch = 16, cex = 0.5)
#  
#  circos.trackPlotRegion(factors = factors, ylim = c(0, 1),bg.col = NA, bg.border = NA)
#  circos.updatePlotRegion(sector.index = "a", bg.border = "black")
#  x1 = runif(100)
#  y1 = runif(100)
#  circos.points(x1, y1, pch = 16, cex = 0.5)
#  
#  circos.trackPlotRegion(factors = factors, ylim = c(0, 1))
#  circos.trackPlotRegion(factors = factors, ylim = c(0, 1))
#  circos.clear()

## ----circlize_part2, echo = FALSE, out.width = "0.8\\textwidth", fig.cap = "Only plot subset of sectors in certain tracks."----
source("src/intro-18-part2.R")

## ----eval=FALSE----------------------------------------------------------
#  par(mar = c(1, 1, 1, 1))
#  factors = letters[1:4]
#  circos.initialize(factors = factors, xlim = c(0, 1))
#  circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
#      circos.text(0.5, 0.5, "outer circos")
#  })
#  circos.clear()
#  
#  par(new = TRUE)
#  circos.par("canvas.xlim" = c(-2, 2), "canvas.ylim" = c(-2, 2))
#  factors = letters[1:3]
#  circos.initialize(factors = factors, xlim = c(0, 1))
#  circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
#      circos.text(0.5, 0.5, "inner circos")
#  })
#  circos.clear()

## ----circlize_nested, echo = FALSE, out.width = "\\textwidth", fig.cap = "An outer circos plot plus an inner one."----
source("src/intro-19-nested.R")

## ----eval=FALSE----------------------------------------------------------
#  par(mar = c(1, 1, 1, 1))
#  factors = letters[1:4]
#  circos.par("canvas.xlim" = c(-1, 1.5), "canvas.ylim" = c(-1, 1.5), start.degree = -45)
#  circos.initialize(factors = factors, xlim = c(0, 1))
#  circos.trackPlotRegion(ylim = c(0, 1), bg.col = NA, bg.border = NA)
#  circos.updatePlotRegion(sector.index = "a")
#  circos.text(0.5, 0.5, "first one")
#  circos.updatePlotRegion(sector.index = "b")
#  circos.text(0.5, 0.5, "first one")
#  
#  circos.clear()
#  
#  par(new = TRUE)
#  circos.par("canvas.xlim" = c(-1.5, 1), "canvas.ylim" = c(-1.5, 1), start.degree = -45)
#  circos.initialize(factors = factors, xlim = c(0, 1))
#  circos.trackPlotRegion(ylim = c(0, 1), bg.col = NA, bg.border = NA)
#  circos.updatePlotRegion(sector.index = "d")
#  circos.text(0.5, 0.5, "second one")
#  circos.updatePlotRegion(sector.index = "c")
#  circos.text(0.5, 0.5, "second one")
#  
#  circos.clear()

## ----circlize_separated, echo = FALSE, out.width = "\\textwidth", fig.cap = "Two separated circos plots"----
source("src/intro-20-seperated.R")

## ----eval=FALSE, echo=TRUE-----------------------------------------------
#  library(circlize)
#  par(mar = c(1, 1, 1, 1))
#  factors = letters[1:4]
#  lim = c(1, 1.1, 1.2, 1.3)
#  for(i in 1:4) {
#      circos.par("canvas.xlim" = c(-lim[i], lim[i]),
#          "canvas.ylim" = c(-lim[i], lim[i]), "track.height" = 0.4)
#      circos.initialize(factors = factors, xlim = c(0, 1))
#      circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA)
#      circos.updatePlotRegion(sector.index = factors[i], bg.border = "black")
#      circos.points(runif(10), runif(10), pch = 16)
#      circos.clear()
#      par(new = TRUE)
#  }
#  par(new = FALSE)

## ----circlize_diff_radius, echo = FALSE, out.width = "0.8\\textwidth", fig.cap = "Sectors with different radius."----
source("src/intro-21-diffradius.R")

## ----eval=FALSE, echo=TRUE-----------------------------------------------
#  set.seed(12345)
#  par(mar = c(1, 1, 1, 1))
#  factors = letters[1:4]
#  circos.par("canvas.xlim" = c(-1.5, 1.5), "canvas.ylim" = c(-1.5, 1.5), "gap.degree" = 10)
#  circos.initialize(factors = factors, xlim = c(0, 1))
#  circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
#      circos.points(1:20/20, 1:20/20)
#  })
#  circos.lines(c(1/20, 0.5), c(1/20, 3), sector.index = "d", straight = TRUE)
#  circos.text(0.5, 3, "mark", sector.index = "d", adj = c(0.5, 0))
#  
#  circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
#      circos.points(1:20/20, 1:20/20)
#  })
#  text(0, 0, "this is\nthe center", cex = 1.5)
#  legend("bottomleft", pch = 1, legend = "this is the legend")
#  circos.clear()

## ----circlize_outside, echo = FALSE, out.width = "0.8\\textwidth", fig.cap = "Draw outside the cell and combine with canvas coordinate."----
source("src/intro-22-outside.R")

## ----circlize_layout, out.width = "\\textwidth", fig.cap = "Arrange multiple circos plots."----
layout(matrix(1:9, 3, 3))
for(i in 1:9) {
    factors = 1:8
    par(mar = c(0.5, 0.5, 0.5, 0.5))
    circos.par(cell.padding = c(0, 0, 0, 0))
    circos.initialize(factors, xlim = c(0, 1))
    circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.05,
        bg.col = rand_color(8), bg.border = NA)
    for(i in 1:20) {
        se = sample(1:8, 2)
        circos.link(se[1], runif(2), se[2], runif(2), 
            col = rand_color(1, transparency = 0.4), border = NA)
    }
    circos.clear()
}

