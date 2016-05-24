### R code from vignette source 'ex-qt-mapping-model.Rnw'

###################################################
### code chunk number 1: qt-mvc-mapper-map
###################################################
data(Cars93, package="MASS")
model <- qdataFrameModel(Cars93, editable=names(Cars93))
mapper <- Qt$QDataWidgetMapper()
mapper$setModel(model)
##
label <- Qt$QLineEdit()
mapper$addMapping(label, 1)


###################################################
### code chunk number 2: qt-mvc-mapper-select
###################################################
table_view <- Qt$QTableView()
table_view$setModel(model)
qconnect(table_view$selectionModel(), "currentRowChanged", 
         function(cur,prev) mapper$setCurrentIndex(cur$row()))


###################################################
### code chunk number 3: qt-mvc-mapper-layout
###################################################
window <- Qt$QWidget()
layout <- Qt$QVBoxLayout()
window$setLayout(layout)
layout$addWidget(table_view)
layout$addWidget(label)


window$show()
