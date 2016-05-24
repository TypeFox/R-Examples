### R code from vignette source 'intro.Rnw'

###################################################
### code chunk number 1: options
###################################################
options(width=72)


###################################################
### code chunk number 2: conventional-scatterplot (eval = FALSE)
###################################################
## x <- rnorm(100)
## plot(x)


###################################################
### code chunk number 3: event-driven-scatterplot (eval = FALSE)
###################################################
## x <- rnorm(100)
## drawScatterplot <- function() { plot(x) }
## # <code that registers drawScatterplot with event loop>
## # control is then passed to event loop


###################################################
### code chunk number 4: drawing-basics-paintFun
###################################################
df <- data.frame(x = rnorm(100, 300, 100),
                 y = rnorm(100, 200, 80))
scatterPainter <- function(layer, painter) {
  qdrawCircle(painter, df$x, df$y, 5)
}


###################################################
### code chunk number 5: drawing-basics-qlayer (eval = FALSE)
###################################################
## library(qtpaint)
## library(qtbase)
## scene <- qscene()
## scatterLayer <- qlayer(scene, paintFun = scatterPainter)


###################################################
### code chunk number 6: drawing-basics-scene-view (eval = FALSE)
###################################################
## view <- qplotView(scene = scene)
## print(view)


###################################################
### code chunk number 7: drawing-basics-glyph (eval = FALSE)
###################################################
## circle <- qglyphCircle()


###################################################
### code chunk number 8: drawing-basics-glyph-paintFun
###################################################
scatterPainter <- function(layer, painter) {
  qdrawGlyph(painter, circle, df$x, df$y)
}


###################################################
### code chunk number 9: drawing-coords-limits (eval = FALSE)
###################################################
## scene <- qscene()
## scatterLayer <- qlayer(scene, paintFun = scatterPainter,
##                        limits = qrect(0, 0, 600, 400))
## print(qplotView(scene = scene))


###################################################
### code chunk number 10: layout-layers (eval = FALSE)
###################################################
## scene <- qscene()
## figLayer <- qlayer(scene)
## titleLayer <- qlayer(figLayer, titlePainter, rowSpan = 2)
## yaxis <- qlayer(figLayer, yAxisPainter, row = 1)
## plotLayer <- qlayer(figLayer, plotPainter, row = 1, col = 1)
## xaxis <- qlayer(figLayer, xAxisPainter, row = 2, col = 1)


###################################################
### code chunk number 11: layout-layers-addLayer (eval = FALSE)
###################################################
## figLayer[0, 0, rowSpan=2] <- titleLayer
## figLayer[1, 0] <- yaxis
## figLayer[1, 1] <- plotLayer
## figLayer[2, 1] <- xaxis


###################################################
### code chunk number 12: layout-preferred-dims (eval = FALSE)
###################################################
## layout <- figLayer$gridLayout()
## layout$setRowPreferredHeight(0, 75)
## layout$setRowPreferredHeight(1, 400)
## layout$setRowPreferredHeight(2, 75)
## layout$setColumnPreferredWidth(0, 75)
## layout$setColumnPreferredWidth(1, 600)


###################################################
### code chunk number 13: layout-stretch-factor (eval = FALSE)
###################################################
## layout$setRowStretchFactor(0, 0)
## layout$setRowStretchFactor(2, 0)
## layout$setColumnStretchFactor(0, 0)


###################################################
### code chunk number 14: input-basics-hoverMoveEvent
###################################################
pointAdder <- function(layer, event) {
  df <<- rbind(df, event$pos())
  qupdate(scene)
}


###################################################
### code chunk number 15: input-basic-qlayer (eval = FALSE)
###################################################
## pointAddingLayer <- qlayer(paintFun = scatterPainter,
##                            mouseReleaseFun = pointAdder)


###################################################
### code chunk number 16: input-mapping-qprimitives
###################################################
pointBrusher <- function(layer, event) {
  rect <- qrect(0, 0, 20, 20)
  mat <- layer$deviceTransform(event)$inverted()
  rect <- mat$mapRect(rect) # 20x20 square now in data space
  pos <- event$pos()
  rect$moveCenter(pos) # centered on the pointer data pos
  hits <- layer$locate(rect) # get indices in rectangle
  df$color[hits] <- "blue" # color points blue
  qupdate(scene)
}


