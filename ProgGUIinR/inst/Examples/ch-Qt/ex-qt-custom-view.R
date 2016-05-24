### R code from vignette source 'ex-qt-custom-view.Rnw'

###################################################
### code chunk number 1: ex-qt-custom-view.Rnw:17-18
###################################################
library(qtbase)


###################################################
### code chunk number 2: CustomView
###################################################
qsetClass("MeanLabel", Qt$QLabel, 
          function(model, column = 0, parent = NULL) 
          {
            super(parent)
            this$model <- model
            this$column <- column
            updateMean()     # initialize text
            qconnect(model, "dataChanged", 
                     function(top_left, bottom_right) {
                       if (top_left$column() <= column && 
                           bottom_right$column() >= column)
                         updateMean()
                     })
          })


###################################################
### code chunk number 3: label
###################################################
qsetMethod("updateMean", MeanLabel, function() {
  if(is.null(model)) {
    text <- "No model"
  } else {
    DF <- qdataFrame(model)
    colname <- colnames(DF)[column + 1L]
    text <- sprintf("Mean for '%s': %s", colname, 
                    mean(DF[,colname]))
  }
  this$text <- text
}, access="private")


###################################################
### code chunk number 4: testItOut
###################################################
model <- qdataFrameModel(mtcars, editable = colnames(mtcars))

table_view <- Qt$QTableView()
table_view$setModel(model)
table_view$setEditTriggers(Qt$QAbstractItemView$DoubleClicked)
##
mean_label <- MeanLabel(model)
##
window <- Qt$QWidget()
layout <- Qt$QVBoxLayout()
window$setLayout(layout)
layout$addWidget(table_view)
layout$addWidget(mean_label)


###################################################
### code chunk number 5: ex-qt-custom-view.Rnw:74-76
###################################################
window$show()
window$raise()


