### =========================================================================
### QtPainter objects
### -------------------------------------------------------------------------
###
if(FALSE){
setRefClass("QtPainter",
            fields = list(
              transform = function(value) {
                if (missing(value))
                  qdeviceTransform(painter)
                else qdeviceTransform(painter) <- value
              },
              strokeColor = function(value) {
                if (missing(value))
                  qstrokeColor(painter)
                else qstrokeColor(painter) <- value
              },
              fillColor = function(value) {
                if (missing(value))
                  qfillColor(painter)
                else qfillColor(painter) <- value
              },
              font = function(value) {
                if (missing(value))
                  qfont(painter)
                else qfont(painter) <- value
              },
              lineWidth = function(value) {
                if (missing(value))
                  qlineWidth(painter)
                else qlineWidth(painter) <- value
              },
              dash = function(value) {
                if (missing(value))
                  qdash(painter)
                else qdash(painter) <- value
              },
              glyphExpansion = function(value) {
                if (missing(value))
                  qglyphExpansion(painter)
                else qglyphExpansion(painter) <- value
              },
              antialias = function(value) {
                if (missing(value))
                  qantialias(painter)
                else qantialias(painter) <- value
              },
              fontMetrics = function(value) {
                if (missing(value))
                  qfontMetrics(painter)
                else qfontMetrics(painter) <- value
              }
              ),
            methods = list(
              drawLine = function(x, y, stroke = NULL) {
                qdrawLine(x, y, stroke)
              },
              drawSegment = function(x0, y0, x1, y1, stroke = NULL) {
                notImplemented("drawSegment")
              },
              drawPoint = function(x, y, stroke = NULL) {
                notImplemented("drawPoint")
              },
              drawRect = function(xleft, ybottom, xright, ytop,
                stroke = NULL, fill = NULL)
              {
                notImplemented("drawRect")
              },
              drawCircle = function(x, y, r, stroke = NULL, fill = NULL) {
                notImplemented("drawCircle")
              },
              drawPolygon = function(x, y, stroke = NULL, fill = NULL) {
                notImplemented("drawPolygon")
              },
              newPath = function() {
                notImplemented("newPath")
              },
              drawPath = function(path, stroke = NULL, fill = NULL) {
                notImplemented("drawPath")
              },
              drawText = function(text, x, y,
                halign = c("center", "left", "right"),
                valign = c("center", "basecenter", "baseline", "bottom", "top"),
                rot = 0, color = NULL, cex = 1.0, hcex = cex, vcex = cex)
              {
                notImplemented("drawText")
              },
              drawImage = function(image, x, y) {
                notImplemented("drawImage")
              },
              drawGlyph = function(glyph, x, y, cex = NULL, stroke = NULL,
                fill = NULL)
              {
                notImplemented("drawGlyph")
              },
              textExtents = function(text) {
                notImplemented("textExtents")
              },
              strWidth = function(text) {
                notImplemented("strWidth")
              },
              strHeight = function(text) {
                notImplemented("strHeight")
              }
              ),
            contains = "Painter")
}
