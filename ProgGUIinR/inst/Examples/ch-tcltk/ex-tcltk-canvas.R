### R code from vignette source 'ex-tcltk-canvas.Rnw'

###################################################
### code chunk number 1: ex-tcltk-canvas.Rnw:8-14
###################################################
## Canvas example of moving a point. See tkcanvas for much more
library(tcltk)
window <- tktoplevel()
tkwm.title(window, "Move canvas object example")
canvas <- tkcanvas(window, width = 450, height = 300)
tkpack(canvas, expand = TRUE, fill = "both")


###################################################
### code chunk number 2: ex-tcltk-canvas.Rnw:20-25
###################################################
x <- 200; y <- 150; r <- 6
item <- tkcreate(canvas, "oval", x - r, y - r, x + r, y + r,
                 width = 1, outline = "black",
                 fill = "blue")
tkaddtag(canvas, "point", "withtag", item)


###################################################
### code chunk number 3: ex-tcltk-canvas.Rnw:32-36
###################################################
tkitembind(canvas, "point", "<Any-Enter>", function()
           tkitemconfigure(canvas, "current", fill = "red"))
tkitembind(canvas, "point", "<Any-Leave>", function()
           tkitemconfigure(canvas, "current", fill = "blue"))


###################################################
### code chunk number 4: ex-tcltk-canvas.Rnw:42-49
###################################################
last_pos <- numeric(2)            # global to track position
tag_selected <- function(W, x, y) {
  tkaddtag(W,  "selected",  "withtag",  "current")
  tkitemraise(W, "current")
  last_pos <<- as.numeric(c(x, y))
}
tkitembind(canvas, "point", "<Button-1>",  tag_selected)


###################################################
### code chunk number 5: moveSelectedPoint
###################################################
move_selected <- function(W, x, y) {
  pos <- as.numeric(c(x,y))
  tkmove(W, "selected", pos[1] - last_pos[1], 
                        pos[2] - last_pos[2])
  last_pos <<- pos
}
tkbind(canvas, "<B1-Motion>", move_selected)


###################################################
### code chunk number 6: ex-tcltk-canvas.Rnw:70-72
###################################################
tkbind(canvas, "<ButtonRelease-1>", 
       function(W) tkdtag(W, "selected"))


