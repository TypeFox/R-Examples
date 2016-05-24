
###################################################
### code chunk number 134: FormLayoutExample
###################################################
window <- Qt$QWidget()
window$setWindowTitle("Wrapper for 'dnorm' function")
window$setLayout(layout <- Qt$QFormLayout())
sapply(c("quantile", "mean", "sd"), function(statistic) {
  layout$addRow(statistic, Qt$QLineEdit())
})
layout$addRow(Qt$QCheckBox("log"))


###################################################
### code chunk number 135: Layouts.Rnw:492-493
###################################################
window$show(); window$raise()
