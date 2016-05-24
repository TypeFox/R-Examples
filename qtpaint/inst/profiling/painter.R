## Performance testing

library(qtpaint)

ops <- c("glyph", "polygon", "rectangle", "point", "line", "segment", "circle")

test <- function(op, opengl, antialias, filled) {
  op <- match(op, ops)
  wid <- .Call("qanviz_painter_test", op - 1L, as.logical(opengl),
               as.logical(antialias), as.logical(filled))
}

params <- expand.grid(antialias = c(TRUE, FALSE), opengl = c(TRUE, FALSE),
                      filled = c(TRUE, FALSE), op = ops)
shapes <- c("glyph", "polygon", "rectangle", "circle")
params <- subset(params, op %in% shapes | (!op %in% shapes & !filled))

## output stored as csv
writeLines("op,ogl,aa,filled,time", "painter.csv")

for (i in seq(nrow(params))) do.call(test, params[i,])

timings <- read.csv("painter-ogl2.csv")
