### R code from vignette source 'ex-qt-edit-data-frame.Rnw'

###################################################
### code chunk number 1: ex-qt-edit-data-frame.Rnw:20-28
###################################################
## Create a model for displaying a data frame -- like qdataFrameModel -- only
## with more Qt-like control over the display of columnms
## The key methods are data(), setData() and flags()
##
## For performance reasons -- which are considerable here and the main reason
## we discourage this approach -- 
## we avoid header Data and put variable names in first row!
library(qtbase)


###################################################
### code chunk number 2: DfModel
###################################################
qsetClass("DfModel", Qt$QAbstractTableModel, 
          function(DF = data.frame(V1 = character(0)), 
                   parent = NULL) 
          {
            super(parent)
            this$DF <- DF
          })


###################################################
### code chunk number 3: ex-qt-edit-data-frame.Rnw:46-50
###################################################
qsetProperty("DF", DfModel, write = function(DF) {
  this$.DF <- DF
  dataChanged(index(0, 0), index(nrow(DF), ncol(DF)))
})


###################################################
### code chunk number 4: ex-qt-edit-data-frame.Rnw:56-60
###################################################
qsetMethod("rowCount", DfModel, 
           function(index) nrow(this$DF) + 1)
qsetMethod("columnCount", DfModel, 
           function(index) ncol(this$DF))


###################################################
### code chunk number 5: displayRoleMethod
###################################################
displayRole <- function(x, row, ...) UseMethod("displayRole")
displayRole.default <- function(x, row) 
  sprintf("%s", x[row])
displayRole.numeric <- function(x, row) 
  sprintf("%.2f", x[row])
displayRole.integer <- function(x, row) 
  sprintf("%d", x[row])


###################################################
### code chunk number 6: ex-qt-edit-data-frame.Rnw:86-125
###################################################
## Not shown
## text alignment to indicate to user different types of numeric values
textAlignmentRole <- function(x, row, context) UseMethod("textAlignmentRole")
textAlignmentRole.default <- function(x, row, context) Qt$Qt$AlignCenter
textAlignmentRole.integer <- function(x, row, context) Qt$Qt$AlignRight
textAlignmentRole.numeric <- function(x, row, context) Qt$Qt$AlignRight

## sets the background color
## returns a QBrush object
backgroundRole <- function(x, row, context) UseMethod("backgroundRole")
backgroundRole.default <- function(x, row, context) Qt$QBrush(Qt$QColor(0,0,0,0)) # transparent
backgroundRole.factor <- function(x, row, context) Qt$QBrush(Qt$QColor("yellow"))

##' Size hint role
##' XXX Doesn't get called
sizeHintRole <- function(x, row, context) UseMethod("sizeHintRole")
sizeHintRole.default <- function(x, row, context) {
  sz <- max(sapply(x, function(i) nchar(as.character(i))))
  avg <- Qt$QFontMetric(Qt$QFont())$averageCharWidth()
  Qt$QSize(sz*avg, sample(c(25L, 50L), 1))
}

## some text for a tooltip
toolTipRole <- function(x, row, context) UseMethod("toolTipRole")
toolTipRole.default <- function(x, row, context)
  sprintf("Value in vector with class %s", paste(class(x), collapse=","))
toolTipRole.factor <- function(x, row, context) {
  x <- levels(x)
  n <- length(x)
  p <- 6
  q <- n %/% p
  r <- n %% p
  l <- paste(paste(apply(matrix(x[1:(p*q)], byrow=T, ncol=6), 1, paste, collapse=","), collapse=",\n"),
        paste(x[(n-r+1):n], collapse=","), sep=",\n")
  sprintf("Factor with levels:\n%s", l)
}
toolTipRole.logical <- function(x, row, context) sprintf("Logical vector")

whatsThisRole <- toolTipRole


###################################################
### code chunk number 7: ex-qt-edit-data-frame.Rnw:130-148 (eval = FALSE)
###################################################
## qsetMethod("data", DfModel, function(index, role) {
##   row <- index$row()
##   col <- index$column() + 1
## 
##   if(role == Qt$Qt$DisplayRole) {
##     if(row > 0)
##       displayRole(DF[,col], row)
##     else
##       names(DF)[col]
##   } else if(role == Qt$Qt$EditRole) {
##     if(row > 0)
##       as.character(DF[row, col])
##     else
##       names(DF)[col]
##   } else {
##     NULL
##   }
## })


###################################################
### code chunk number 8: data
###################################################
## this is not shown in text, but is the definition of the data method
qsetMethod("data", DfModel, function(index, role) {
  row <- index$row()
  col <- index$column() + 1

  
  if(role == Qt$Qt$DisplayRole) {
    if(row > 0)
      displayRole(DF[,col], row)
    else
      names(DF)[col]
  } else if(role == Qt$Qt$EditRole) {
    if(row > 0)
      as.character(DF[row, col])
    else
      names(DF)[col]
  } else if(role == Qt$Qt$TextAlignmentRole) {
    if(row > 0)
      textAlignmentRole(DF[, col], row)
    else
      Qt$Qt$AlignCenter | Qt$Qt$AlignTop
  } else if(role == Qt$Qt$BackgroundRole) {
    if(row > 0)
      backgroundRole(DF[, col], row)
    else
      Qt$QBrush(Qt$QColor("gray"))
  } else if(role == Qt$Qt$ForegroundRole) {
    if(row == 0)
      Qt$QBrush(Qt$QColor("white"))
  } else if(role == Qt$Qt$ToolTipRole) {
    if(row > 0)
      toolTipRole(DF[,col], row)
    else
      ""
  } else if(role == Qt$Qt$WhatsThisRole) {
    if(row > 0)
      whatsThisRole(DF[,col], row)
    else
      ""
  } else if(role == Qt$Qt$SizeHintRole) {
    if(row > 0)
      sizeHintRole(DF[,col], row)
    else
      NULL
  } else {
    NULL
  }
})


###################################################
### code chunk number 9: ex-qt-edit-data-frame.Rnw:205-213
###################################################
qsetMethod("flags", DfModel, function(index) {
  if(!index$isValid()) {
    return(Qt$Qt$ItemIsEnabled)
  } else {
    current_flags <- super("flags", index)
    return(current_flags | Qt$Qt$ItemIsEditable)
  }
})


###################################################
### code chunk number 10: fitIn
###################################################
fitIn <- function(x, value) UseMethod("fitIn")
fitIn.default <- function(x, value) value
fitIn.numeric <- function(x, value) as.numeric(value)


###################################################
### code chunk number 11: notShown
###################################################
## more methods for fit in
fitIn.integer <- function(x, value) as.integer(value)
fitIn.logical <- function(x, value) {
  if(toupper(value) %in% c("T","F","TRUE","FALSE")) {
    as.logical(value)
  } else {
    as.logical(as.numeric(value))
  }
}


###################################################
### code chunk number 12: ex-qt-edit-data-frame.Rnw:241-261
###################################################
qsetMethod("setData", DfModel, function(index, value, role) {
  if(index$isValid() && role == Qt$Qt$EditRole) {
    DF <- this$DF
    row <- index$row()
    col <- index$column() + 1

    if(row > 0) {
      x <- DF[, col]
      DF[row, col] <- fitIn(x, value)
    } else {
      names(DF)[col] <- value
    }
    this$DF <- DF
    dataChanged(index, index)

    return(TRUE)
  } else {
     super("setData", index, value, role)
  }
})


###################################################
### code chunk number 13: ex-qt-edit-data-frame.Rnw:267-280
###################################################
qsetMethod("setColumn", DfModel, function(col, value) {
  ## pad with NA if needed
  n <- nrow(this$DF)
  if(length(value) < n)
    value <- c(value, rep(NA, n - length(value)))
  value <- value[1:n]
  DF <- this$DF
  DF[,col] <- value
  this$DF <- DF         # only notify about this column
  dataChanged(index(0, col - 1), 
              index(rowCount() - 1, col - 1))
  return(TRUE)
})


###################################################
### code chunk number 14: addColumn
###################################################
qsetMethod("addColumn", DfModel, function(name, value) {
  DF <- this$DF
  if(name %in% names(DF)) {
    return(setColumn(min(which(name == names(DF))), value))
  }  
  beginInsertColumns(Qt$QModelIndex(),
                     columnCount(), columnCount())
  DF[[name]] <- value
  this$DF <- DF
  endInsertColumns()
  return(TRUE)
})


###################################################
### code chunk number 15: ex-qt-edit-data-frame.Rnw:302-306
###################################################
model <- DfModel(mtcars)

view <- Qt$QTableView()
view$setModel(model)


###################################################
### code chunk number 16: customizeView
###################################################
trigger_flag <- Qt$QAbstractItemView$DoubleClicked | 
                Qt$QAbstractItemView$SelectedClicked |
                Qt$QAbstractItemView$EditKeyPressed
view$setEditTriggers(trigger_flag)
view$verticalHeader()$setHidden(TRUE)
view$horizontalHeader()$setHidden(TRUE)


###################################################
### code chunk number 17: ex-qt-edit-data-frame.Rnw:319-321
###################################################
view$show()
view$raise()


