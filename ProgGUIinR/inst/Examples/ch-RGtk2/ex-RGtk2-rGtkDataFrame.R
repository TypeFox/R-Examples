### R code from vignette source 'ex-RGtk2-rGtkDataFrame.Rnw'

###################################################
### code chunk number 1: ex-RGtk2-rGtkDataFrame.Rnw:2-3
###################################################
library(RGtk2)


###################################################
### code chunk number 2: callBackEdit
###################################################
editCallBack <- function(cell, path, arg3, ...) {
  if(nargs() == 3) {
    userData <- arg3; newValue <- NA    # no newValue (toggle)
  } else {
    newValue <- arg3; userData = ..1    # ..1 is first component of ...
  }
  rGtkStore <- userData$view$getModel() 
  i <- as.numeric(path) + 1
  j <- userData$column
  newValue <- try(switch(userData$type,
                         "integer" = as.integer(as.numeric(newValue)),
                         "character" = as.character(newValue),
                         "numeric" = as.numeric(newValue),
                         "factor"  = as.character(newValue),
                         "logical" =  !as.logical(rGtkStore[i,j])),
                  silent=TRUE)
  
  if(inherits(newValue,"try-error")) {
    sprintf("Failed to coerce new value to type %s",userData$type)
    return(FALSE)
  }
  
  if(userData$type == "factor") {
    curLevels <- levels(rGtkStore[,j])
    if(! newValue %in% curLevels) {
      cat(gettext("Can't add level to a factor."))
      return(FALSE)
    }
  }
  
  rGtkStore[i,j] <- newValue            # assign value
  return(FALSE)
}


###################################################
### code chunk number 3: AddColumnWithType
###################################################
gtkTreeViewAddColumnWithType <-
  function(view,
           name="",
           type=c("rowname","numeric","integer","character",
             "logical","factor","icon"),
           viewCol,                     # 1-based column of view
           storeCol                     # 1-based column for rGtkDataFrame
           ) {
    
    type = match.arg(type)
    
    ## define the cell renderer
    cr <- switch(type,
                 "logical" = gtkCellRendererToggle(),
                 "factor" = gtkCellRendererCombo(),
                 gtkCellRendererText())
    
    ## the new column we will add
    vc <- gtkTreeViewColumn()
    vc$packStart(cr, TRUE)
    vc$setTitle(name)
    vc$setResizable(TRUE); vc$setClickable(TRUE)
    view$InsertColumn(vc, viewCol - 1)  # viewCol is 1-based
    
    ## add attributes
    switch(type,
           "logical" =  vc$addAttribute(cr, "active",storeCol - 1),
           vc$addAttribute(cr, "text",storeCol - 1)
           )
    if(type == "numeric") cr['xalign'] <- 1
    
    ## set editable/activatable property
    switch(type,
           "logical" = cr["activatable"] <- TRUE,
           cr["editable"] <- TRUE)
    
    if(type == "factor") {              # combo box needs a data store
      cstore <- gtkListStore("gchararray")
      rGtkstore <- view$getModel()
      vals <- rGtkstore[,storeCol, drop=TRUE]
      for(i in as.character(unique(vals))) {
        iter <- cstore$append()
        cstore$setValue(iter$iter,column=0, i)
      }
      cr['model'] <- cstore
      cr['text-column'] <- 0        
    }

    
    ## connect callback to edited/toggled signal
    QT <- gSignalConnect(cr, signal =
                         if(type != "logical") "edited" else "toggled",
                         f = editCallBack, 
                         data = list(view=view,type=type,column=storeCol))
  }


###################################################
### code chunk number 4: keyNav
###################################################
### -- bug with this when not editing
gtkTreeViewAddKeyNavigations <- function(view) {
  ## keyMotionHandler example.
  connectSignal(view,"key-release-event",
                f = function(view, event, userData,...) {
                  
                  keyval = event$getKeyval()
                  cursor = view$getCursor()
                  ## i,j are current positions,
                  i = as.numeric(cursor$path$toString()) + 1
                  vc = cursor[['focus.column']] ## might be focus_column!!
                  j = which(sapply(view$getColumns(), function(i) i == vc))
                  d = dim(view$getModel()) # rGtkStore method
                  
                  setCursorAtCell <- function(view, i, j) {
                    path <- gtkTreePathNewFromString(i-1)
                    vc <- view$getColumn(j - 1)
                    view$setCursor(path=path, focus.column=vc, start.editing=TRUE)
                  }
                  
                  if(keyval == GDK_Return) {
                    ## what do do with return?
                  } else if(keyval == GDK_Up) {
                    setCursorAtCell(view,max(1, i - 1), j)
                  } else if(keyval == GDK_Down) {
                    if(i < d[1]) 
                      setCursorAtCell(view,i + 1, j)
                  } else if(keyval == GDK_Tab) {
                    if(j < d[2]) 
                      setCursorAtCell(view,i, j + 1)
                  }
                },
                data=list(view = view)
                )

}


###################################################
### code chunk number 5: testIt
###################################################
DF = data.frame(
  logical = c(TRUE, TRUE, FALSE),
  character = c("one","two","three"),
  factor = factor(c("ctrl","trt1","trt2")),
  integer = 1:3,
  numeric = rnorm(3),
  stringsAsFactors=FALSE)


###################################################
### code chunk number 6: ex-RGtk2-rGtkDataFrame.Rnw:166-168
###################################################
store <- rGtkDataFrame(DF)
view <- gtkTreeView(store)


###################################################
### code chunk number 7: ex-RGtk2-rGtkDataFrame.Rnw:172-178
###################################################
nms <- names(DF)
QT <- sapply(1:ncol(DF), function(i) {
  type <- class(DF[,i])[1]
  view$addColumnWithType(name = nms[i], type, viewCol = i, storeCol = i)
  
})


###################################################
### code chunk number 8: ex-RGtk2-rGtkDataFrame.Rnw:182-184
###################################################
vc <- gtkTreeViewColumn()
newColNo <- view$insertColumn(vc, -1)           # -1 = end


###################################################
### code chunk number 9: AddNavigations
###################################################
ID <- view$addKeyNavigations()


###################################################
### code chunk number 10: PackWidget
###################################################
sw <- gtkScrolledWindow()
sw$setPolicy("automatic","automatic")
sw$add(view)

w <- gtkWindow(); w$setTitle("rGtkDataFrame example")
w$add(sw)


