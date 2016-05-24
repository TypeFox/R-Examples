textlabel <- function(x, y = NULL, labels = seq_along(x), bg.col = grey(0.85), adj = NULL, pos = NULL, offset = 0.5, vfont = NULL, cex = 1, col = NULL, font = NULL, ...){
##
## Plot text with rectangular background, i.e. text labels
## NOTE: Not all arguments above have been dealt with....
#
## Width and height of text 
.width <- strwidth(labels, cex = cex, font = font, vfont = vfont, ...) * 1.4
.height <- strheight(labels, cex = cex, font = font, vfont = vfont, ...) * 1.8
#
## Corners of rectangle
.coords <- xy.coords(x, y)
.coords.rec <- c(xleft = .coords$x - .width/2, ybottom = .coords$y - .height, xright = .coords$x + .width/2, ytop = .coords$y + .height)
#
## 
rect(xleft = .coords.rec["xleft"], ybottom = .coords.rec["ybottom"], xright = .coords.rec["xright"], ytop = .coords.rec["ytop"], col = bg.col, border = "black")
text(x, y, labels = labels, cex = cex, font = font)
#
##
return(invisible(.coords.rec))
}

