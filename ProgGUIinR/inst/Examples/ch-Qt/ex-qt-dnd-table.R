### R code from vignette source 'ex-qt-dnd-table.Rnw'

###################################################
### code chunk number 1: ex-qt-dnd-table.Rnw:1-3
###################################################
library(qtbase)
#rm(list=ls())                           # clear out


###################################################
### code chunk number 2: qt-mvc-dnd-variable-selector
###################################################
qsetClass("VariableSelector", Qt$QWidget, 
          function(parent = NULL) {
  super(parent)
  ## widgets
  this$df_combo_box <- Qt$QComboBox()
  this$variable_list <- Qt$QListView()
  this$variable_list$setModel(
       qdataFrameModel(data.frame(), this, 
                       useRoles = TRUE))
  this$variable_list$dragEnabled <- TRUE

  ## layout
  layout <- Qt$QVBoxLayout()
  layout$addWidget(df_combo_box)
  layout$addWidget(variable_list)
  variable_list$setSizePolicy(Qt$QSizePolicy$Expanding, 
                        Qt$QSizePolicy$Expanding)
  setLayout(layout)
  
  updateDataSets()
  qconnect(df_combo_box, "activated(int)", function(ind) {
    this$dataFrame <- df_combo_box$currentText
  })
})


###################################################
### code chunk number 3: qt-mvc-dnd-update-datasets
###################################################
qsetMethod("updateDataSets", VariableSelector, function() {
  current_text <- df_combo_box$currentText
  df_combo_box$clear()
  DFs <- ProgGUIinR:::avail_dfs(.GlobalEnv)
  if(length(DFs)) {
    this$df_combo_box$addItems(DFs)
    if(is.null(current_text) || !current_text %in% DFs) {
      this$df_combo_box$currentIndex <- -1
      this$dataFrame <- NULL
    } else {
      this$df_combo_box$currentIndex <- 
        which(current_text == DFs)
      this$dataFrame <- current_text
    }
  }
})


###################################################
### code chunk number 4: ex-qt-dnd-table.Rnw:72-83
###################################################
## The \function{getIcon} helper function provides an icon from the class of a
## column:
require(grid)
getIcon <- function(x) {
  f <- tempfile()
  png(file=f, width=16, height=16)
  grid::grid.newpage()
  grid::grid.draw(ProgGUIinR:::make_icon(x))
  dev.off()
  Qt$QIcon(f)
}


###################################################
### code chunk number 5: qt-mvc-dnd-dataset
###################################################
qsetProperty("dataFrame", VariableSelector, 
             write = function(DF) {
               if (is.null(DF))
                 DF <- data.frame()
               else if (is.character(DF)) 
                 DF <- get(DF, .GlobalEnv)
               ##
               model <- variable_list$model()
               icons <- lapply(DF, getIcon)
               qdataFrame(model) <- 
                 data.frame(variable=names(DF),
                            variable.decoration=I(icons))
               this$.dataFrame <- DF
               dataFrameChanged()
             })


###################################################
### code chunk number 6: qt-mvc-dnd-datasetChanged
###################################################
qsetSignal("dataFrameChanged", VariableSelector)


###################################################
### code chunk number 7: DropLabelRotation
###################################################
qsetClass("VariableLabel", Qt$QLabel, function(parent=NULL) {
  super(parent)
  this$rotation <- 0L
  setAcceptDrops(TRUE)
  setAlignment(Qt$Qt$AlignHCenter | Qt$Qt$AlignVCenter)
})


###################################################
### code chunk number 8: qt-mvc-dnd-rotation
###################################################
qsetProperty("rotation", VariableLabel)
qsetProperty("variable_name", VariableLabel)


###################################################
### code chunk number 9: qt-mvc-dnd-drop
###################################################
qsetSignal("variableNameDropped", VariableLabel)


###################################################
### code chunk number 10: qt-mvc-dnd-get-variable-name
###################################################
variableNameFromMimeData <- function(mime_data) {
  name <- NULL
  RDA_MIME_TYPE <- "application/x-rlang-transport"
  if(mime_data$hasFormat(RDA_MIME_TYPE)) {
    name_list <- unserialize(mime_data$data(RDA_MIME_TYPE))
    if (length(name_list) && is.character(name_list[[1]]))
      name <- name_list[[1]]
  }
  name
}


###################################################
### code chunk number 11: ex-qt-dnd-table.Rnw:159-170
###################################################
qsetMethod("dragEnterEvent", VariableLabel, function(event) {
  mime_data <- event$mimeData()
  if(!is.null(variableNameFromMimeData(mime_data))) {
    setForegroundRole(Qt$QPalette$Dark)
    event$acceptProposedAction()
  }
})
qsetMethod("dragLeaveEvent", VariableLabel, function(event) {
  setForegroundRole(Qt$QPalette$WindowText)
  event$accept()
})


###################################################
### code chunk number 12: dropEvent
###################################################
qsetMethod("dropEvent", VariableLabel, function(event) {
  setForegroundRole(Qt$QPalette$WindowText)  
  mime_data <- event$mimeData()
  this$variable_name <- variableNameFromMimeData(mime_data)
  if(!is.null(variable_name)) {
    this$text <- variable_name
    variableNameDropped()
    setBackgroundRole(Qt$QPalette$Window)
    event$acceptProposedAction()
  }
})


###################################################
### code chunk number 13: ex-qt-dnd-table.Rnw:195-208
###################################################
qsetMethod("paintEvent", VariableLabel, function(event) {
  painter <- Qt$QPainter()
  painter$begin(this)
  
  painter$save()
  painter$translate(width / 2, height / 2)
  painter$rotate(-(rotation))
  rect <- painter$boundingRect(0, 0, 0, 0, 
                               Qt$Qt$AlignCenter, text)
  painter$drawText(rect, Qt$Qt$AlignCenter, text)
  painter$restore()
  painter$end()
})


###################################################
### code chunk number 14: XtabsWidget
###################################################
qsetClass("XtabsWidget", Qt$QWidget, function(parent = NULL) {
  super(parent)
  initWidgets()
  initLayout()
})


###################################################
### code chunk number 15: initWidgets
###################################################
qsetMethod("initWidgets", XtabsWidget, function() {
  this$xlabel <- VariableLabel()
  qconnect(xlabel, "variableNameDropped", invokeXtabs)

  this$ylabel <- VariableLabel()
  pt <- ylabel$font$pointSize()
  ylabel$minimumWidth <- 2*pt; ylabel$maximumWidth <- 2*pt
  ylabel$rotation <- 90L
  qconnect(ylabel, "variableNameDropped", invokeXtabs)
  
  this$table_view <- Qt$QTableView()
  table_view$setModel(qdataFrameModel(data.frame(), this))
  clearLabels()
})


###################################################
### code chunk number 16: ex-qt-dnd-table.Rnw:253-262
###################################################
## Not shown
qsetMethod("clearLabels", XtabsWidget, function() {
  point_size <- xlabel$font$pointSize()
  xlabel$text <- "Drop x factor here"
  xlabel$minimumHeight <- 2 * point_size
  
  ylabel$text <- "Drop y factor here"
  ylabel$minimumWidth <- 2 * point_size
})


###################################################
### code chunk number 17: ex-qt-dnd-table.Rnw:265-276
###################################################
## Not shown
qsetMethod("initLayout", XtabsWidget, function() {
  layout <- Qt$QGridLayout()
  setLayout(layout)
  layout$addWidget(xlabel, 0, 1, 1, 3)
  layout$addWidget(ylabel, 1, 0, 1, 1)
  layout$addWidget(table_view, 1, 1, 1, 3)
  
  layout$setColumnStretch(2, 1)
  layout$setRowStretch(1, 1)
})


###################################################
### code chunk number 18: ex-qt-dnd-table.Rnw:281-297
###################################################
## Hide call to xtabs
## Return NULL if not okay, otherwise a table object
call_xtabs <- function(DF, x, y) {
  if(is.character(DF))
    DF <- get(DF)
  if(is.null(x)) {
    table <- NULL
  } else if(is.null(y)) {
    f <- formula(sprintf("~ %s", x))
    table <- xtabs(f, data = DF)
  } else { 
    f <- formula(sprintf("~ %s + %s", y, x))
    table <- xtabs(f, data = DF)
  } 
  table
}


###################################################
### code chunk number 19: makeTable
###################################################
qsetMethod("invokeXtabs", XtabsWidget, function() {
  if (is.null(dataFrame))
    return()

  x <- xlabel$variable_name
  y <- ylabel$variable_name
  
  if(!is.null(table <- call_xtabs(dataFrame, x, y)))
     updateTableView(table)
})


###################################################
### code chunk number 20: updateTableWidget
###################################################
qsetMethod("updateTableView", XtabsWidget, function(table) {
  model <- table_view$model()
  if (length(dim(table)) == 1)
    qdataFrame(model) <- data.frame(count = unclass(table))
  else qdataFrame(model) <- data.frame(unclass(table))
})


###################################################
### code chunk number 21: qt-mvc-dnd-dataframe-xtabs
###################################################
qsetProperty("dataFrame", XtabsWidget, 
             write = function(dataFrame) { 
               clearLabels()
               this$.dataFrame <- dataFrame
             })


###################################################
### code chunk number 22: ex-qt-dnd-table.Rnw:336-338
###################################################
## Not shown
require(MASS); data(Cars93); data(Aids2)


###################################################
### code chunk number 23: ex-qt-dnd-table.Rnw:340-349
###################################################
w <- Qt$QSplitter()
w$setWindowTitle("GUI for xtabs()")
w$addWidget(vs <- VariableSelector())
w$addWidget(tw <- XtabsWidget())
w$setStretchFactor(1, 1)
qconnect(vs, "dataFrameChanged", function() {
  tw$dataFrame <- vs$dataFrame
})
w$show(); w$raise()


