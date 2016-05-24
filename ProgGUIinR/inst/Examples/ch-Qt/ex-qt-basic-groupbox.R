
###################################################
### code chunk number 137: groupBoxExample
###################################################
## Not shown, Group box example of title, alignment, flat, checkable
f <- Qt$QGroupBox("Example group box")
lyt <- Qt$QVBoxLayout(); f$setLayout(lyt)
lyt$addWidget(changeTitle <- Qt$QPushButton("Change title"))
lyt$addWidget(changeAlignment <- Qt$QPushButton("Cycle Alignment"))
lyt$addWidget(toggleFlat <- Qt$QPushButton("Toggle flat"))
lyt$addWidget(toggleCheckable <- Qt$QPushButton("Toggle checkable"))
f$show()

qconnect(changeTitle, "clicked", function(checked) {
  f$setTitle("New title")
})
qconnect(changeAlignment, "clicked", function(checked) {
  aligns <- c(Qt$Qt$AlignLeft, Qt$Qt$AlignHCenter, Qt$Qt$AlignRight)
  curAlign <- f$alignment
  ind <-   which(curAlign == aligns)
  f$setAlignment(aligns[c(2,3,1)[ind]])
})
qconnect(toggleFlat, "clicked", function(checked) {
  f$setFlat(!f$flat)
})
qconnect(toggleCheckable, "clicked", function(checked) {
  f$setCheckable(!f$checkable)
})

