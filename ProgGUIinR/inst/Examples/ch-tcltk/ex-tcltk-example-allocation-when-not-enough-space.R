###################################################
### code chunk number 74: NotEnoughSpace
###################################################
## example of last in first covered
## Create this GUI, then shrink window with the mouse
w <- tktoplevel()
f <- ttkframe(w); tkpack(f, expand = TRUE, fill = "both")

g1 <- ttkframe(f); tkpack(g1, expand = TRUE, fill = "both")
g2 <- ttkframe(f); tkpack(g2, expand = TRUE, fill = "both")

b11 <- ttkbutton(g1, text = "first")
b12 <- ttkbutton(g1, text = "second")
b21 <- ttkbutton(g2, text = "first")
b22 <- ttkbutton(g2, text = "second")

tkpack(b11, side = "left"); tkpack(b12, side = "left")
tkpack(b21, side="right"); tkpack(b22, side = "right")
## Now shrink window with the mouse, what keeps the size?
