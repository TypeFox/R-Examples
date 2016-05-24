### R code from vignette source 'ex-RGtk2-select-variables.Rnw'

###################################################
### code chunk number 1: VariableSelectionExample
###################################################
## Example showing implementation of variable selection widget where two tables show possible selections
## and selection. Similar to SPSS widget
## Illustrates filtered models, icons in view column
library(RGtk2)


###################################################
### code chunk number 2: ex-RGtk2-select-variables.Rnw:24-25
###################################################
DF <- get(data(Cars93, package="MASS"))


###################################################
### code chunk number 3: ex-RGtk2-select-variables.Rnw:41-43
###################################################
library(ProgGUIinR)                     # for make_icon


###################################################
### code chunk number 4: make_icon
###################################################
make_icon_pixmap <- function(x, ...) {
  require(grid); require(cairoDevice)
  pixmap <- gdkPixmap(drawable = NULL, width = 16, height=16, 
                      depth = 24)
  asCairoDevice(pixmap)
  grid.newpage()
  grid.draw(ProgGUIinR:::make_icon(x))
  dev.off()
  gdkPixbufGetFromDrawable(NULL, pixmap, NULL, 0,0,0,0,-1,-1)
}


###################################################
### code chunk number 5: model
###################################################
model_df <- data.frame(Variables = I(sort(names(DF))),
                       icon = I(sapply(DF, make_icon_pixmap)),
                       selected = rep(FALSE, ncol(DF)))
model <- rGtkDataFrame(model_df)


###################################################
### code chunk number 6: filterModels
###################################################
selected_filter <- model$filter()
selected_filter$setVisibleColumn(2)
unselected_filter <- model$filter()
unselected_filter$setVisibleFunc(function(model, iter) {
  !model$get(iter, 2)[[1]]
})


###################################################
### code chunk number 7: views
###################################################
views <- list()
views$unselected_view <- gtkTreeView(unselected_filter)
views$selected_view <- gtkTreeView(selected_filter)
##
sapply(views, function(view) {
  selection <- view$getSelection()
  selection$setMode('multiple')
})


###################################################
### code chunk number 8: viewColumns
###################################################
make_view_column <- function() {
  column <- gtkTreeViewColumn()
  column$setTitle("Variable")
  column$packStart(cell_renderer <- gtkCellRendererPixbuf())
  column$addAttribute(cell_renderer, "pixbuf", 1L)
  column$packStart(cell_renderer <- gtkCellRendererText())
  column$addAttribute(cell_renderer, "text", 0L)
  column
}
sapply(views, function(view) 
       view$insertColumn(make_view_column(), 0))


###################################################
### code chunk number 9: extendAPI
###################################################
## add to the gtkTreeView API for convenience
gtkTreeViewSelectedIndices <- function(object) {
  model <- object$getModel()          # Filtered!
  paths <- object$getSelection()$getSelectedRows()$retval
  path_strings <- sapply(paths, function(i) {
    model$convertPathToChildPath(i)$toString()
  })
  if(length(path_strings) == 0)
    integer(0)
  else
    as.numeric(path_strings) + 1 # 1-based
}
## does object have selection?
gtkTreeViewHasSelection <-
  function(obj) length(obj$selectedIndices()) > 0


###################################################
### code chunk number 10: buttons
###################################################
buttons <- list()
buttons$unselect_button <- gtkButton("<")
buttons$select_button <- gtkButton(">")
toggleSelectionOnClick <- function(button, view) {
  gSignalConnect(button, "clicked", function(button) {
    message("clicked")
    ind <- view$selectedIndices()
    model[ind, "selected"] <- !model[ind, "selected"]
  })
}
sapply(1:2, function(i) toggleSelectionOnClick(buttons[[i]], 
                                               views[[3-i]]))


###################################################
### code chunk number 11: sensitiveButtons
###################################################
sapply(buttons, gtkWidgetSetSensitive, FALSE)
trackSelection <- function(button, view) {
  gSignalConnect(view$getSelection(), "changed", 
     function(x) button['sensitive'] <- view$hasSelection())
}
sapply(1:2, function(i) trackSelection(buttons[[i]], 
                                       views[[3-i]]))


###################################################
### code chunk number 12: guiLayout
###################################################
window <- gtkWindow(show=FALSE)
window$setTitle("Select variables example")
window$setDefaultSize(600, 400)
hbox <- gtkHBox()
window$add(hbox)
## scrollwindows
scrolls <- list()
scrolls$unselected_scroll <- gtkScrolledWindow()
scrolls$selected_scroll <- gtkScrolledWindow()
mapply(gtkContainerAdd, object = scrolls, widget = views)
mapply(gtkScrolledWindowSetPolicy, scrolls, 
       "automatic", "automatic")
## buttons
button_box <- gtkVBox()
centered_box <- gtkVBox()
button_box$packStart(centered_box, expand=TRUE, fill = FALSE)
centered_box$setSpacing(12)
sapply(buttons, centered_box$packStart, expand = FALSE)
##
hbox$packStart(scrolls$unselected_scroll, expand = TRUE)
hbox$packStart(button_box, expand = FALSE)
hbox$packStart(scrolls$selected_scroll, expand = TRUE)


###################################################
### code chunk number 13: packButtons
###################################################
window$show()


