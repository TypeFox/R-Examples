
###################################################
### code chunk number 142: notShown
###################################################
library(gWidgets)
options(guiToolkit = "RGtk2")


###################################################
### code chunk number 143: ourData
###################################################
data("Cars93", package = "MASS")


###################################################
### code chunk number 144: layout
###################################################
window <- gwindow("Spotfire example", visible = FALSE)
notebook <- gnotebook(cont = window)


###################################################
### code chunk number 145: page1
###################################################
descr <- glabel(gettext("A basic GUI to explore a data set"), 
                cont = notebook, label = gettext("About"))


###################################################
### code chunk number 146: page2
###################################################
group <- ggroup(cont = notebook, label = gettext("Explore..."))
left_group <- ggroup(cont = group, horizontal = FALSE) 
right_group <- ggroup(cont = group, horizontal = FALSE)


###################################################
### code chunk number 147: ex-gWidgets-spotfire.Rnw:60-61
###################################################
ggraphics(cont = left_group)


###################################################
### code chunk number 148: ex-gWidgets-spotfire.Rnw:72-78
###################################################
tbl <- gtable(Cars93, cont = left_group, multiple = TRUE, 
              filter.FUN = "manual")
size(tbl) <- c(500, 200)                # set size
label_group <- ggroup(cont = left_group)
addSpring(label_group)
no_cases <- glabel("", cont = label_group)


###################################################
### code chunk number 149: filters
###################################################
filter_frame <- gframe(gettext("Filter by:"), 
                       cont = right_group, expand = TRUE)


###################################################
### code chunk number 150: filterLayout
###################################################
lyt <- glayout(cont = filter_frame)
widget_list <- list() # store widgets
lyt[1,1] <- "Type:"
lyt[1,2] <- (widget_list$Type <- 
             gcombobox(c("", levels(Cars93$Type)), 
                       cont = lyt))

lyt[2,1] <- "Cylinders:"
lyt[2,2] <- (widget_list$Cylinders <- 
             gcombobox(c("", levels(Cars93$Cyl)), cont = lyt))


###################################################
### code chunk number 151: handlers
###################################################
update_data_frame <- function(...) {
  vals <- lapply(widget_list, svalue)
  vals <- vals[vals != ""] 
  out <- sapply(names(vals), function(i) {
    Cars93[[i]] == vals[[i]]
  })
  ind <- apply(out, 1, function(x) Reduce("&&", x))
  ## update table
  visible(tbl) <- ind
  ## update label
  nsprintf <- function(n, msg1, msg2,...)
    ngettext(n, sprintf(msg1, n), sprintf(msg2,n), ...)
  svalue(no_cases) <- nsprintf(sum(ind),"%s case", "%s cases")
}


###################################################
### code chunk number 152: ex-gWidgets-spotfire.Rnw:155-162
###################################################
update_graphic <- function(...) {
  ind <- visible(tbl)
  if(any(ind))
    plot(MPG.city ~ Weight, data = Cars93[ind,])
  else
    plot.new()
}


###################################################
### code chunk number 153: applyHandler
###################################################
callback <- function(h, ...) {
  update_data_frame()
  update_graphic()
}
sapply(widget_list, addHandlerChanged, handler = callback)


###################################################
### code chunk number 154: tableHandler
###################################################
addHandlerClicked(tbl, handler = function(h,...) {
  update_graphic()
  ind <- svalue(h$obj, index = TRUE)
  points(MPG.city ~ Weight, cex = 2, col = "red", pch = 16, 
         data = Cars93[ind,])
})


###################################################
### code chunk number 155: initialGraphic
###################################################
visible(window) <- TRUE
update_graphic()
