## convenient constructor wrappers for some simple types
## 

qrect <- function(x0, y0, x1, y1) {
  if (missing(x0)) # a 'null' rectangle
    return(Qt$QRectF())  
  if (length(x0) == 4L) {
    y0 <- x0[,2]
    x0 <- x0[,1]
  }
  if (length(x0) == 2L) {
    x1 <- x0[2]
    x0 <- x0[1]
  }
  if (length(y0) == 2L) {
    y1 <-  y0[2]
    y0 <-  y0[1]
  }
  if (is.integer(c(x0, y0, x1, y1)))
    cons <- Qt$QRect
  else cons <- Qt$QRectF
  cons(x0, y0, (x1 - x0), (y1 - y0))
}

qpoint <- function(x, y) {
  if (missing(x)) # 'null' point
    return(Qt$QPointF())
  if (length(x) == 2) {
    y <- x[2]
    x <- x[1]
  }
  if (is.integer(c(x, y)))
    cons <- Qt$QPoint
  else cons <- Qt$QPointF
  cons(x, y)
}

qsize <- function(width, height) {
  if (missing(width)) # 'null' size
    return(Qt$QSizeF())
  if (length(width) == 2) {
    height <- width[2]
    width <- width[1]
  }
  if (is.integer(c(width, height)))
    cons <- Qt$QSize
  else cons <- Qt$QSizeF
  cons(width, height)
}

qpolygon <- function(x = NULL, y = NULL) {
  xy <- xy.coords(x, y)
  if (is.integer(c(x, y))) {
    pointCon <- Qt$QPoint
    con <- Qt$QPolygon
  } else {
    pointCon <- Qt$QPointF
    con <- Qt$QPolygonF
  }
  con(mapply(pointCon, xy$x, xy$y))
}

qfont <- function(family = baseFont$family(), pointsize = baseFont$pointSize(),
                  weight = baseFont$weight(),
                  italic = baseFont$style() == Qt$QFont$StyleItalic,
                  baseFont = Qt$QApplication$font())
{
  Qt$QFont(family, pointsize, weight, italic)
}

qpen <- function(brush = qbrush(), width = 0L, style = Qt$Qt$SolidLine,
                 cap = Qt$Qt$SquareCap, join = Qt$Qt$BevelJoin)
{
  if (!is(brush, "QBrush"))
    brush <- qbrush(brush)  
  Qt$QPen(brush, width, style, cap, join)
}

qbrush <- function(color = qcolor(), style = Qt$Qt$SolidPattern)
{
  if (!is(color, "QColor"))
    color <- qcolor(color)
  Qt$QBrush(color, style)
}

qcolor <- function(red = 0, green = 0, blue = 0, alpha = 255)
{
    if (is.character(red) || is.matrix(red))
    {
        rgbvals <-
            if (is.character(red)) col2rgb(red, alpha = TRUE)[,1]
        else
            red[,1]
        red <- rgbvals["red"]
        green <- rgbvals["green"]
        blue <- rgbvals["blue"]
        if (missing(alpha)) alpha <- rgbvals["alpha"]
    }
    Qt$QColor(red, green, blue, alpha)
}

qtransform <- function(m11 = 1.0, m12 = 0.0, m13 = 0.0, m21 = 0.0, m22 = 1.0,
                       m23 = 0.0, m31 = 0.0, m32 = 0.0, m33 = 1.0)
{
  if (is.matrix(m11)) {
    stopifnot(is.numeric(m11) && ncol(m11) == 3L && nrow(m11) == 3L)
    m11 <- m11[1]; m12 <- m11[4]; m13 <- m11[7]
    m21 <- m11[2]; m22 <- m11[5]; m23 <- m11[8]
    m31 <- m11[3]; m32 <- m11[6]; m33 <- m11[9]
  }
  Qt$QTransform(m11, m12, m13, m21, m22, m23, m31, m32, m33)
}

as.QImage <- function(x, ...) UseMethod("as.QImage")
as.QImage.default <- function(x, ...) {
  if (!is.matrix(x) || !is.integer(x) || (nrow(x) != 4 && nrow(x) != 3))
    rgb <- col2rgb(x, TRUE)
  else rgb <- x
  Qt$QImage(as.raw(rgb), ncol(x), nrow(x), ncol(x) * nrow(rgb),
            if (nrow(rgb) == 3) Qt$QImage$Format_RGB888
            else Qt$QImage$Format_ARGB32)
}
