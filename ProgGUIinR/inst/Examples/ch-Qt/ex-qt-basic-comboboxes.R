
###################################################
### code chunk number 215: Widgets.Rnw:456-460
###################################################
df <- data.frame(name=state.name, region=state.region,
                 population=state.x77[,'Population'], 
                 stringsAsFactors=FALSE)
states_by_region <- split(df, df$region)


###################################################
### code chunk number 216: QComboBox
###################################################
state_combo <- Qt$QComboBox()
region_combo <- Qt$QComboBox()
region_combo$addItems(names(states_by_region))


###################################################
### code chunk number 217: qt-widget-combo-currentitem
###################################################
region_combo$currentText
region_combo$currentIndex                     # 0-based


###################################################
### code chunk number 218: qt-widget-combo-clear-index
###################################################
region_combo$currentIndex <- -1


###################################################
### code chunk number 219: Widgets.Rnw:491-495
###################################################
qconnect(region_combo, "activated(int)", function(index) {
  state_combo$clear()
  state_combo$addItems(states_by_region[[index+1]]$name)
})


###################################################
### code chunk number 220: Widgets.Rnw:502-509
###################################################
window <- Qt$QGroupBox("Two combo boxes")
layout <- Qt$QFormLayout()
window$setLayout(layout)
layout$addRow("Region:", region_combo)
layout$addRow("State:", state_combo)
layout$fieldGrowthPolicy <-  # grow combo boxes
  Qt$QFormLayout$AllNonFixedFieldsGrow


###################################################
### code chunk number 221: Widgets.Rnw:512-513 (eval = FALSE)
###################################################
window$show(); window$raise()
