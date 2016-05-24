
###################################################
### code chunk number 120: stretch_expand
###################################################
## Not shown
## Example of various combinations of stretch, alignment, sizePolicy
## no stretch, no expand
w <- Qt$QWidget(); w$setWindowTitle("No stretch, no expand")
w$setLayout(g <- Qt$QHBoxLayout())
buttons <- sapply(letters[1:3], Qt$QPushButton)
sapply(buttons, g$addWidget)
w$setMinimumSize(400, 50)
w$show()

## no stretch, expand for first
w <- Qt$QWidget();  w$setWindowTitle("No stretch, no expand, size Policy")
w$setLayout(g <- Qt$QHBoxLayout())

buttons <- sapply(letters[1:3], Qt$QPushButton)
sapply(buttons, g$addWidget)

b <- g$itemAt(0L)$widget()
b$setSizePolicy(Qt$QSizePolicy$Expanding, Qt$QSizePolicy$Fixed) 
for(i in 1:2) {
  b <- g$itemAt(i)$widget()
  b$setSizePolicy(Qt$QSizePolicy$Fixed, Qt$QSizePolicy$Fixed) 
}
w$setMinimumSize(400, 50)
w$show()                

## stretch
w <- Qt$QWidget(); w$setWindowTitle("Using stretch factors")
w$setLayout(g <- Qt$QHBoxLayout())
buttons <- sapply(letters[1:3], Qt$QPushButton)
for(i in 1:3) g$addWidget(buttons[[i]], stretch=i)


w$setMinimumSize(400, 50)
w$show()                

## no stretch; alignment
w <- Qt$QWidget(); w$setWindowTitle("Using alignment")
w$setLayout(g <- Qt$QHBoxLayout())
g$addWidget(Qt$QPushButton("NW"), stretch=0, Qt$Qt$AlignLeft | Qt$Qt$AlignTop)
g$addWidget(Qt$QPushButton("Center"), stretch=0, Qt$Qt$AlignHCenter | Qt$Qt$AlignVCenter)
g$addWidget(Qt$QPushButton("SE"), stretch=0, Qt$Qt$AlignRight | Qt$Qt$AlignBottom)
w$setMinimumSize(400, 400)
w$show()
##
