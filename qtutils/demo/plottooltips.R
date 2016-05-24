library(qtutils)

rscene <- qsceneDevice()
(view <- QT(rscene))
view$setDragMode(Qt$QGraphicsView$NoDrag)

## A convoluted way to add tooltips to a plot

plot(state.x77[, "Illiteracy"], state.x77[, "Income"],
     bty = "n", axes = FALSE, xlab = "", ylab = "")

pitems <- rscene$items() # the points in the plot

for (i in seq_along(pitems))
    pitems[[i]]$setToolTip(rownames(state.x77)[i])

axis(1); axis(2); box()
title(main = "U.S.States",
      xlab = "Percent Illiterate (1970, population)",
      ylab = "Per capita income (1974)")



