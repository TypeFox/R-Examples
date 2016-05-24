### R code from vignette source 'ex-RGtk2-UImanager.Rnw'

###################################################
### code chunk number 1: not-shown
###################################################
## sample RGtk2 menu
library(RGtk2)


###################################################
### code chunk number 2: ex-RGtk2-UImanager.Rnw:15-16
###################################################
uimanager = gtkUIManager()


###################################################
### code chunk number 3: ex-RGtk2-UImanager.Rnw:22-25
###################################################
someAction <- function(action,...) 
  statusbar$push(statusbar$getContextId("message"), 
                 action$getName())


###################################################
### code chunk number 4: ex-RGtk2-UImanager.Rnw:29-30
###################################################
Quit <- function(...) win$destroy()


###################################################
### code chunk number 5: Define-first-action-group
###################################################
firstActionGroup <- gtkActionGroup("firstActionGroup")
firstActionEntries <- list(
  ## name,ID,label,accelerator,tooltip,callback
  file = list("File",NULL,"_File",NULL,NULL,NULL),
  new = list("New", "gtk-new", "_New", "<control>N", 
    "New document", someAction),
  sub = list("Submenu", NULL, "S_ub", NULL, NULL, NULL),
  open = list("Open", "gtk-open", "_Open", "<ctrl>0", 
    "Open document", someAction),
  save = list("Save", "gtk-save", "_Save", "<alt>S", 
    "Save document", someAction),
  quit = list("Quit", "gtk-quit", "_Quit", "<ctrl>Q", 
    "Quit", Quit),
  edit = list("Edit", NULL, "_Edit", NULL, NULL, NULL),
  undo = list("Undo", "gtk-undo", "_Undo", "<ctrl>Z", 
    "Undo change", someAction),
  redo = list("Redo", "gtk-redo", "_Redo", "<ctrl>U", 
    "Redo change", someAction)
)


###################################################
### code chunk number 6: "Insert action group"
###################################################
firstActionGroup$addActions(firstActionEntries)
uimanager$insertActionGroup(firstActionGroup, 0) # 0-based


###################################################
### code chunk number 7: How-to-set-sensitivity
###################################################
redo <- firstActionGroup$getAction("Redo")
redo['sensitive'] <- FALSE


###################################################
### code chunk number 8: ex-RGtk2-UImanager.Rnw:86-93
###################################################
helpActionGroup <- gtkActionGroup("helpActionGroup")
helpActionEntries <- list(
  help = list("Help", "", "_Help", "", "", NULL),
  about = list("About", "gtk-about", "_About", "", "", 
    someAction)
)
helpActionGroup$addActions(helpActionEntries)


###################################################
### code chunk number 9: "a toggle action"
###################################################
toggleActions <- list(
  tooltips = list("UseTooltips", NULL, "Use _Tooltips", "<control>T", 
    "Enable tooltips", someAction, TRUE)
)
helpActionGroup$addToggleActions(toggleActions)


###################################################
### code chunk number 10: "insert help action group"
###################################################
uimanager$insertActionGroup(helpActionGroup, 1)


###################################################
### code chunk number 11: "Load UI from file"
###################################################
xml_file <- system.file("resources", "ex-menus.xml", package="ProgGUIinR")
id <- uimanager$addUiFromFile(xml_file)


###################################################
### code chunk number 12: "Retrieve menubar and toolbar from the uimanager"
###################################################
menubar <- uimanager$getWidget("/menubar")
toolbar <- uimanager$getWidget("/toolbar")


###################################################
### code chunk number 13: "define statusbar"
###################################################
statusbar <- gtkStatusbar()


###################################################
### code chunk number 14: Define-window-add-accelerator-group
###################################################
win <- gtkWindow(show=TRUE)
win$setTitle("Window example")
accelgroup <- uimanager$getAccelGroup()
win$addAccelGroup(accelgroup)


###################################################
### code chunk number 15: setup-box
###################################################
box <- gtkVBox()
win$add(box)

box$packStart(menubar, expand=FALSE, fill=FALSE, 0)
box$packStart(toolbar, expand=FALSE, fill= FALSE, 0)
contentArea <- gtkVBox()
box$packStart(contentArea, expand=TRUE, fill=TRUE, 0)
contentArea$packStart(gtkLabel("Content Area"))
box$packStart(statusbar, expand=FALSE, fill=FALSE, 0)


