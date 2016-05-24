### R code from vignette source 'ex-RGtk2-combobox-entry.Rnw'

###################################################
### code chunk number 1: ex-RGtk2-combobox-entry.Rnw:1-4
###################################################
## a combobox that learns as you go.
## no tooltip per item, but here we add as detail
library(RGtk2)


###################################################
### code chunk number 2: ex-RGtk2-combobox-entry.Rnw:14-18
###################################################
model <- rGtkDataFrame(data.frame(filename = character(0), 
                                  visits = character(0), 
                                  nvisits = integer(0), 
                                  stringsAsFactors = FALSE))


###################################################
### code chunk number 3: ex-RGtk2-combobox-entry.Rnw:32-34
###################################################
combo_box <- gtkComboBoxEntryNewWithModel(model, 
                                          text.column = 0)


###################################################
### code chunk number 4: ConfigureCellRenderers
###################################################
cell_renderer <- gtkCellRendererText()
combo_box$packStart(cell_renderer)
combo_box$addAttribute(cell_renderer, "text", 1)
cell_renderer['foreground'] <- "gray50"
cell_renderer['ellipsize'] <- "end"
cell_renderer['style'] <- "italic"
cell_renderer['alignment'] <- "right"


###################################################
### code chunk number 5: helperFunction2
###################################################
callHelpFunction <- function(combo_box, value) {
  model <- combo_box$getModel()
  ind <- match(value, model[,1,drop=TRUE])
  nvisits <- model[ind, "nvisits"] <- model[ind, "nvisits"]+1
  model[ind, "visits"] <- 
    sprintf(ngettext(nvisits,"%s visit","%s visits"), nvisits)
  ## select for easier editing
  combo_box$getChild()$selectRegion(start = 0, end = -1)
  help(value)
}
gSignalConnect(combo_box, "changed", 
               f = function(combo_box, ...) {
                 if(combo_box$getActive() >= 0) {
                   value <- combo_box$getActiveText()
                   callHelpFunction(combo_box, value)
                 }
               })


###################################################
### code chunk number 6: ex-RGtk2-combobox-entry.Rnw:98-110
###################################################
gSignalConnect(combo_box$getChild(), "activate", 
       f = function(combo_box, entry, ...) {
         value <- entry$getText()
         if(!any(value == combo_box$getModel()[,1])) {
           model <- combo_box$getModel()
           tmp <- data.frame(filename = value, visits = "", 
                             nvisits = 0, 
                             stringsAsFactors = FALSE)
           model$appendRows(tmp)
         }
         callHelpFunction(combo_box, value)
       }, data = combo_box, user.data.first = TRUE)


###################################################
### code chunk number 7: Layout
###################################################
window <- gtkWindow(show = FALSE)
window['border-width'] <- 15
hbox <- gtkHBox(); window$add(hbox)
hbox$packStart(gtkLabel("Help on:"))
hbox$packStart(combo_box, expand = TRUE, fill = TRUE)
#
window$show()


