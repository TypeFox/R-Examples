## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 7)

## ------------------------------------------------------------------------
library(orsifronts)
cols <- hcl(seq(0, 240, length = nrow(orsifronts)), c = 50)
plot(orsifronts, col = cols, lwd = 2)
degAxis(1)
degAxis(2)
box()

## ------------------------------------------------------------------------
as.data.frame(orsifronts)
plot(orsifronts, xlim = c(60, 180), col = cols, asp = 1/cos(55 * pi / 180), lwd = 4)
legend("topleft", sprintf("%s (%s)", orsifronts$front, orsifronts$name),  
       col = cols, lwd = 4, cex = 0.8, bty = "n")
degAxis(1)
degAxis(2)
box()

