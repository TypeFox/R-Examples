plotColorScale <- function(
    ##title<< Add a color scale to plots
    ##description<< plotColorScale is a wrapper function around color.legend to ease its usage.
   col               ##<< vector of color strings defining the palette to use
   , zlim = c()      ##<< numeric vector (of length 2) defining the upper and lower limit
                     ##   of the values mapped to the color scale.
   , pos = list(x = c(1.02,1.08), y = c(0.1,0.9))
                     ##<< list of x and y coordinates (relative) defining the lower left and
                     ##   upper right edge of the scale.
   , align = 'rb'    ##<< character: alignment option passed to color.legend
   , gradient = 'y'  ##<< character: orientation option passed to color.legend
   , cex = 1         ##<< numeric: character expansion factor for the text labels
   , cex.title = 1
   , title = 'cts/px'      ##<< character: the title of the color scale
   , outer.range = c(FALSE, FALSE) ##<< logical: whether to extend the scale over the zlim borders
                     ##   at its bottom and top.
  , legend =  c()
)
 ##seealso<<
  ##\code{\link{color.legend}}

{
 if (class(pos) == 'character') {
   if (pos == 'right') {
     pos = list(x=c(1.01,1.03), y=c(0.1,0.9))
     gradient = 'y'
   } else if (pos == 'bottom') {  
     pos = list(x=c(0.1,0.9), y=c(-0.15,-0.1))
     gradient = 'x'
   } else if (pos == 'bottom_in') {  
     pos = list(x=c(0.1,0.9), y=c(0.1,0.15))
     gradient = 'x'
   } else if (pos == 'top_in') {  
     pos = list(x=c(0.1,0.9), y=c(0.9,0.95))
     gradient = 'x'
   }
 }  
 coords <-       userCoords(pos$x, pos$y)
 if (sum(outer.range) > 0) {
   pos.outer = pos
   if (gradient == 'y') {
     coords$y <- pos.outer$y - 0.05 * outer.range * c(-1, 1)
   } else {
     coords$x <- pos.outer$x - 0.05 * outer.range * c(-1, 1)
   }
   legend.ext <- paste(c('<', '>'), zlim)
   legend.ext[!outer.range] <- '  '
   coords.outer <-       userCoords(pos.outer$x, pos.outer$y)

   color.legend(coords.outer$x[1],coords.outer$y[1],coords.outer$x[2],coords.outer$y[2],
                rect.col= rep(colorChangeDarkness(col[c(1, length(col))], c(0.5, 0.5)), each = 10),
                legend = legend.ext, gradient = gradient, align = align, cex = cex)
 }
  if (length(zlim) == 0 & length(legend) == 0) {
    legend = ' '
  } else if (length(legend) == 0) {
    legend = seq(zlim[1],zlim[2], length.out = 5)
  }
 color.legend(coords$x[1],coords$y[1],coords$x[2],coords$y[2], rect.col = col,
              legend = legend, gradient=gradient,align=align,cex=cex)
 if (nchar(title) > 0) {
   pos.text = userCoords(mean(pos$x), y = max(pos$y))
   if (gradient == 'y') {
     if( outer.range[2]) {
       pos.text = userCoords(mean(pos$x), y = max(pos$y) + 0.05)
     } else {
            pos.text = userCoords(max(pos$x), y = max(pos$y))
     }
   }
   par(xpd=NA)
   text(labels = title, x = pos.text$x, y = pos.text$y, pos = 3, cex = cex.title)
 }
 ##value<< Nothing is returned.
} 
