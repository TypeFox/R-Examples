### R code from vignette source 'ex-RGtk2-options-in-table.Rnw'

###################################################
### code chunk number 1: editableTableForCollectingOptions
###################################################
## GUI for configuring options -- in a table
library(RGtk2)


###################################################
### code chunk number 2: ex-RGtk2-options-in-table.Rnw:21-28
###################################################
opts <- c("main", "sub", "xlab", "ylab", "line", "outer")
DF <- data.frame(option = opts,
           value = c("", "", "", "", "0", "FALSE"),
           class = c(rep("character",4),"integer", "logical"),
           edit_color = rep("gray95", 6),
           dirty = rep(FALSE, 6),
           stringsAsFactors = FALSE)


###################################################
### code chunk number 3: model
###################################################
model <- rGtkDataFrame(DF)
view <- gtkTreeView(model)
##
cell_renderer <- gtkCellRendererText()
cell_renderer['background'] <- 'gray80'
view$insertColumnWithAttributes(position = -1,
                                title = "Option",
                                cell = cell_renderer,
                                text = 1 - 1)


###################################################
### code chunk number 4: secondColumn
###################################################
cell_renderer <- gtkCellRendererText()
cell_renderer['editable'] <- TRUE
view$insertColumnWithAttributes(position = -1,
                                title = "Value",
                                cell = cell_renderer,
                                text = 2 - 1,
                                background = 4 - 1
                                )


###################################################
### code chunk number 5: editConnect
###################################################
gSignalConnect(cell_renderer, "edited", 
    function(cell_renderer, path, new.text, user.data) {
      model <- user.data$model
      i <- as.numeric(path) + 1; j <- user.data$column
      val <- as(new.text, model[i, 'class'])
      model[i,j] <- as(val, "character")   
      model[i, 'dirty'] <- TRUE                 # mark dirty
      model[i, 'edit_color'] <- 'gray70'        # change color
    }, data=list(model=model, column=2))


###################################################
### code chunk number 6: ex-RGtk2-options-in-table.Rnw:85-92
###################################################
window <- gtkWindow(show=FALSE)
window['title'] <- "Option editor"
window$setSizeRequest(300,500)
scrolled_window <- gtkScrolledWindow()
window$add(scrolled_window)
scrolled_window$add(view)
window$show()


###################################################
### code chunk number 7: ex-RGtk2-options-in-table.Rnw:114-120
###################################################
require(helpr, quietly=TRUE)
package <- "graphics"; topic <- "title"
rd <- helpr:::parse_help(helpr:::pkg_topic(package, topic), 
                         package = package)
descs <- rd$params$args
names(descs) <- sapply(descs, function(i) i$param)


###################################################
### code chunk number 8: ex-RGtk2-options-in-table.Rnw:129-143
###################################################
view["has-tooltip"] <- TRUE
gSignalConnect(view, "query-tooltip", 
       function(view, x, y, key_mode, tooltip, user.data) {
         out <- view$getTooltipContext(x, y, key_mode)
         if(out$retval) {
           model <- view$getModel()
           i <- as.numeric(out$path$toString()) + 1
           val <- model[i, "option"]
           txt <- descs[[val]]$desc
           txt <- gsub("code>","b>", txt)  # no code in Pango
           tooltip$setMarkup(txt)
         }
         out$retval
       })


