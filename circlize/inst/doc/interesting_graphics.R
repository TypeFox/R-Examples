## ----echo = FALSE--------------------------------------------------------
library(knitr)
opts_chunk$set(fig.pos = "", fig.align = "center")

library(circlize)
circos.initialize = function(...) {
    circos.par(unit.circle.segments = 300)
    circlize::circos.initialize(...)
}

## ----interesting_clock, fig.align = "center", out.width = "0.6\\textwidth", fig.cap = "A clock."----
library(circlize)
factors = "a"  # any name is OK
circos.par(gap.degree = 0, cell.padding = c(0, 0, 0, 0), start.degree = 90)
circos.initialize(factors = factors, xlim = c(0, 12))
circos.trackPlotRegion(factors = factors, ylim = c(0, 1), bg.border = NA)
circos.axis(sector.index = "a", major.at = 0:12, labels = "",
    direction = "inside", major.tick.percentage = 0.3)
circos.text(1:12, rep(0.5, 12), 1:12, facing = "downward")
arrows(0, 0, 0, 0.7)    
arrows(0, 0, 0.4, 0)

circos.clear()

## ----interesting_dartboard, fig.align = "center", out.width = "0.6\\textwidth", fig.cap = "A dartboard."----
factors = 1:20  # just indicate there are 20 sectors
circos.par(gap.degree = 0, cell.padding = c(0, 0, 0, 0),
    start.degree = 360/20/2, track.margin = c(0, 0), clock.wise = FALSE)
circos.initialize(factors = factors, xlim = c(0, 1))
circos.trackPlotRegion(ylim = c(0, 1), factors = factors, bg.col = "black",
    track.height = 0.15)
circos.trackText(rep(0.5, 20), rep(0.5, 20),
    labels = c(13, 4, 18, 1, 20, 5, 12, 9, 14, 11, 8, 16, 7, 19, 3, 17, 2, 15, 10, 6),
    factors = factors, col = "#EEEEEE", font = 2, facing = "downward")
circos.trackPlotRegion(ylim = c(0, 1), factors = factors,
    bg.col = rep(c("#E41A1C", "#4DAF4A"), 10), bg.border = "#EEEEEE", track.height = 0.05)
circos.trackPlotRegion(ylim = c(0, 1), factors = factors,
    bg.col = rep(c("black", "white"), 10), bg.border = "#EEEEEE", track.height = 0.275)
circos.trackPlotRegion(ylim = c(0, 1), factors = factors,
    bg.col = rep(c("#E41A1C", "#4DAF4A"), 10), bg.border = "#EEEEEE", track.height = 0.05)
circos.trackPlotRegion(ylim = c(0, 1), factors = factors, 
    bg.col = rep(c("black", "white"), 10), bg.border = "#EEEEEE", track.height = 0.375)
draw.sector(center = c(0, 0), start.degree = 0, end.degree = 360,
    rou1 = 0.1, col = "#4DAF4A", border = "#EEEEEE")
draw.sector(center = c(0, 0), start.degree = 0, end.degree = 360,
    rou1 = 0.05, col = "#E41A1C", border = "#EEEEEE")

circos.clear()

## ----interesting_bagua, fig.align = "center", out.width = "0.6\\textwidth", fig.cap = "Bagua and Taiji."----
factors = 1:8
circos.par(start.degree = 22.5, gap.degree = 6)
circos.initialize(factors = factors, xlim = c(0, 1))

# yang yao is -
add_yang_yao = function() {
    circos.rect(0,0,1,1, col = "black")
}

# yin yao is --
add_yin_yao = function() {
    circos.rect(0,0,0.45,1, col = "black")
    circos.rect(0.55,0,1,1, col = "black")
}
circos.trackPlotRegion(ylim = c(0, 1), factors = factors, bg.border = NA,
    panel.fun = function(x, y) {
        i = get.cell.meta.data("sector.numeric.index")
        if(i %in% c(2, 5, 7, 8)) add_yang_yao() else add_yin_yao()
}, track.height = 0.1)

circos.trackPlotRegion(ylim = c(0, 1), factors = factors, bg.border = NA,
    panel.fun = function(x, y) {
        i = get.cell.meta.data("sector.numeric.index")
        if(i %in% c(1, 6, 7, 8)) add_yang_yao() else add_yin_yao()
    }, track.height = 0.1)

circos.trackPlotRegion(ylim = c(0, 1), factors = factors, bg.border = NA, 
    panel.fun = function(x, y) {
        i = get.cell.meta.data("sector.numeric.index")
        if(i %in% c(4, 5, 6, 7)) add_yang_yao() else add_yin_yao()
    }, track.height = 0.1)

# the bottom of the most recent track
r = get.cell.meta.data("cell.bottom.radius") - 0.1
# draw taiji, note default order is clock wise for `draw.sector`
draw.sector(center = c(0, 0), start.degree = 90, end.degree = -90,
    rou1 = r, col = "black", border = "black")
draw.sector(center = c(0, 0), start.degree = 270, end.degree = 90,
    rou1 = r, col = "white", border = "black")
draw.sector(center = c(0, r/2), start.degree = 0, end.degree = 360,
    rou1 = r/2, col = "white", border = "white")
draw.sector(center = c(0, -r/2), start.degree = 0, end.degree = 360,
    rou1 = r/2, col = "black", border = "black")
draw.sector(center = c(0, r/2), start.degree = 0, end.degree = 360,
    rou1 = r/8, col = "black", border = "black")
draw.sector(center = c(0, -r/2), start.degree = 0, end.degree = 360,
    rou1 = r/8, col = "white", border = "white")

circos.clear()

