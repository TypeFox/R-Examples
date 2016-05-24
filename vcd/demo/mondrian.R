library(vcd)

## shape
foo1 <- c(3, 7, 3, 1.5)
foo2 <- c(2, 6.5, 1.5)
foo <- outer(foo1/sum(foo1), foo2/sum(foo2), "*")

## color
mondrian <- rep("#EAE6E3", 12)
mondrian[1] <- "#DE1024"
mondrian[3] <- "#FFD83B"
mondrian[12] <- "#032349"

## plot
## best visualized with resized display, e.g. using:
## get(getOption("device"))(width = 4.9, height = 7.5)
grid.newpage()
grid.rect(gp = gpar(fill = 1))

mondrianMosaic <- function(x, fill)
  mosaic(x, gp = gpar(col = rep(0, length(fill)), fill = fill),
         legend = FALSE, margins = 0,  newpage = FALSE, keep_aspect_ratio = FALSE)

mondrianMosaic(foo, mondrian)

