## tests/demonstrations for the interactive canvas

library(qtpaint)

options(warn=2)
options(error=recover)

circle <- qglyphCircle()

##circle <- qglyphCircle(r = 0.04)
paths <- lapply(seq(0, by=0.1, length=10), circle$translated, 0)

n <- 100000
x <- rnorm(n) #, 50, 25)
y <- rnorm(n) #, 50, 25)
df <- data.frame(X = x, Y = y)

##data(mtcars)
##df <- mtcars[,c("mpg", "hp")]

##data(iris)
##df <- iris

##fill <- col2rgb(rgb(1, seq(0, 1, length=nrow(df)), 0, 0.5), TRUE)
##fill <- col2rgb(rgb(1, 0, 0, 0.5), TRUE)
fill <- rep("black", nrow(df))
fill <- col2rgb(fill, alpha=TRUE)
##fill[5] <- "blue"
scatterplot <- function(item, painter) {
  ##qstrokeColor(painter) <- NA
  ##qfillColor(painter) <- fill
  ##qantialias(painter) <- FALSE
  ##qdrawText(painter, "x", df[,1], df[,2], color = fill)
  ##qdrawPoint(painter, df[,1], df[,2], stroke = fill)
  ##qdrawCircle(painter, df[,1], df[,2], 5, fill = fill)
  ##qdrawPath(painter, paths, fill = fill)
  qdrawGlyph(painter, circle, df[,1], df[,2], fill=fill)
}

labeled <- rep(FALSE, nrow(df))
labeler <- function(item, painter) {
  mat <- qdeviceTransform(painter)$inverted()
  off <- qvmap(mat, c(5, 5)) - qvmap(mat, c(0, 0))
  df <- df[labeled,]
  qdrawText(painter, rownames(df), df[,1]+off[1], df[,2]+off[2], "left",
            "bottom")
}

margin <- 5
adjust <- c(margin, -margin)

adjustPoint <- Qt$QPointF(margin, margin)

axes <- function(item, painter) {
  qfont(painter) <- qfont(pointsize=12)
  pos <- as.matrix(item$geometry) + adjust
  qdrawText(painter, colnames(df)[1], pos[2], pos[4], "right", "bottom")
  qdrawText(painter, colnames(df)[2], pos[1], pos[3], "left", "top")
}

## pointAdder <- function(item, event) {
##   df <<- rbind(df, event$pos())
##   qupdate(scene)
## }

pointIdentifier <- function(item, event) {
  off <- 20
  rect <- qrect(0, 0, off*2, off*2)
  mat <- item$deviceTransform(event)$inverted()
  rect <- mat$mapRect(rect)
  pos <- event$pos()
  rect$moveCenter(pos)
  hits <- item$locate(rect) + 1L
  hitmat <- as.matrix(df[hits,])
  posmat <- matrix(pos, ncol=2)
  labeled <<- rep(FALSE, nrow(df))
  labeled[hits][Biobase::matchpt(posmat, hitmat)[,1]] <<- TRUE
  qupdate(labels)
}

## boundsPainter <- function(item, painter) {
##   lims <- dim(item)
##   qstrokeColor(painter) <- "red"
##   qdrawRect(painter, lims[1,1], lims[1,2], lims[2,1], lims[2,2])
## }

scene <- qscene()
##qbackground(scene) <- "black"
root <- qlayer(scene)
lims <- qrect(range(df[,1]), range(df[,2]))
points <- qlayer(root, scatterplot,
                 hoverMove = pointIdentifier,
                 cache = TRUE,
                 limits = lims)
labels <- qlayer(root, labeler, cache = FALSE, limits = lims)
##bounds <- qlayer(NULL, boundsPainter)
##qaddGraphicsWidget(root, bounds, 1, 1)
view <- qplotView(scene = scene, opengl = TRUE)
overlay <- view$overlay()
axesOverlay <- qlayer(overlay, axes, cache = TRUE)
print(view)
