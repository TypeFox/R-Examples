## An example of an application window (chapter 16)
## not all the code is "hooked" up.

require(qtbase)
require(qtutils)

###################################################
### code chunk number 422: qt-app-construct
###################################################
main_window <- Qt$QMainWindow()


###################################################
### code chunk number 423: qt-app-table
###################################################
data(mtcars)
model <- qdataFrameModel(mtcars, editable = TRUE)
table_view <- Qt$QTableView()
table_view$setModel(model)
main_window$setCentralWidget(table_view)


###################################################
### code chunk number 424: qt-app-action-construct
###################################################
open_action <- Qt$QAction("Open", main_window)


###################################################
### code chunk number 425: qt-app-action-status-tip
###################################################
open_action$statusTip <- "Load a spreadsheet from a CSV file"


###################################################
### code chunk number 426: qt-app-action-icon
###################################################
style <- Qt$QApplication$style()
open_action$icon <- 
  style$standardIcon(Qt$QStyle$SP_DialogOpenButton)


###################################################
### code chunk number 427: qt-app-action-trigger
###################################################
qconnect(open_action, "triggered", function() {
  filename <- Qt$QFileDialog$getOpenFileName()
  table_view$model <- 
    qdataFrameModel(read.csv(filename), editable = TRUE)
})


###################################################
### code chunk number 428: qt-app-action-toggle
###################################################
save_on_exit_action <- Qt$QAction("Save on exit", main_window)
save_on_exit_action$checkable <- TRUE


###################################################
### code chunk number 429: qt-app-action-radio
###################################################
just_group <- Qt$QActionGroup(main_window)
just_action <- list()
just_action$left <- Qt$QAction("Left Align", just_group)
just_action$right <- Qt$QAction("Right Align", just_group)
just_action$center <- Qt$QAction("Center", just_group)
sapply(just_action, function(action) action$checkable <- TRUE)


###################################################
### code chunk number 430: ApplicationWindow.Rnw:130-137
###################################################
sapply(just_action, function(action) 
       qconnect(action, "changed", function() {
         button_number <- 
           which(sapply(just_action, `[[`, "checked"))
         message("Button ", button_number, " was depressed")
       })
       )


###################################################
### code chunk number 431: actionGroupSignal
###################################################
qconnect(just_group, "triggered", function(action) {
  message(action$text)
})


###################################################
### code chunk number 432: qt-app-action-shortcut
###################################################
open_action$setShortcut(Qt$QKeySequence(Qt$QKeySequence$Open))


###################################################
### code chunk number 433: qt-app-menubar
###################################################
menubar <- Qt$QMenuBar()
main_window$setMenuBar(menubar)


###################################################
### code chunk number 434: qt-app-menubar-addMenu
###################################################
file_menu <- Qt$QMenu("File")
menubar$addMenu(file_menu)
edit_menu <- Qt$QMenu("Edit")
menubar$addMenu(edit_menu)


###################################################
### code chunk number 435: qt-app-menubar-populate
###################################################
file_menu$addAction(open_action)
file_menu$addSeparator()
file_menu$addAction(save_on_exit_action)
file_menu$addSeparator()
quit_action <- file_menu$addAction("Quit")

just_menu <- edit_menu$addMenu("Justification")
sapply(just_action, just_menu$addAction)


###################################################
### code chunk number 436: qt-app-context-add-widgets
###################################################
sort_menu <- Qt$QMenu("Sort by")
sapply(colnames(qdataFrame(model)), sort_menu$addAction)
table_view$addAction(sort_menu$menuAction())


###################################################
### code chunk number 437: qt-app-context-policy
###################################################
table_view$contextMenuPolicy <- Qt$Qt$ActionsContextMenu


###################################################
### code chunk number 438: qt-app-context-signal
###################################################
showCompletionPopup <- function(event, edit) {
  popup <- Qt$QMenu()
  completions <- utils:::matchAvailableTopics(ed$text)
  completions <- head(completions, 10) # trim if large
  sapply(completions, function(completion) {
    action <- popup$addAction(completion)
    qconnect(action, "triggered", 
             function(...) edit$setText(completion))
  })
  popup$popup(edit$mapToGlobal(qpoint(0L,0L)))
}
##
edit <- Qt$QLineEdit()
edit$contextMenuPolicy <- Qt$Qt$CustomContextMenu
qconnect(edit, "customContextMenuRequested", 
         showCompletionPopup, user.data = edit)


###################################################
### code chunk number 439: qt-app-toolbar
###################################################
toolbar <- Qt$QToolBar()
main_window$addToolBar(toolbar)


###################################################
### code chunk number 440: ApplicationWindow.Rnw:336-343
###################################################
## Before adding some actions to our toolbar, we define a function
## \code{getIcon} that loads a \class{QIcon} from a file in the
## \pkg{gWidgets} package:
getIcon <- function(nm) {
  fname <- system.file(sprintf("images/%s.gif", nm), package="gWidgets")
  Qt$QIcon(fname)
}


###################################################
### code chunk number 441: ApplicationWindow.Rnw:348-359
###################################################
file_actions <- list()
file_actions$open <- Qt$QAction("Open", main_window)
file_actions$open$setIcon(getIcon("open"))
file_actions$save <- Qt$QAction("Save", main_window)
file_actions$save$setIcon(getIcon("save"))

plot_actions <- list()
plot_actions$barplot <- Qt$QAction("Barplot", main_window)
plot_actions$barplot$setIcon(getIcon("barplot"))
plot_actions$boxplot <- Qt$QAction("Boxplot", main_window)
plot_actions$boxplot$setIcon(getIcon("boxplot"))


###################################################
### code chunk number 442: ApplicationWindow.Rnw:364-367
###################################################
sapply(file_actions, toolbar$addAction)
toolbar$addSeparator()
sapply(plot_actions, toolbar$addAction)


###################################################
### code chunk number 443: ApplicationWindow.Rnw:369-371
###################################################
main_window$show()
main_window$raise()


###################################################
### code chunk number 444: qt-app-toolbar-style
###################################################
toolbar$setToolButtonStyle(Qt$Qt$ToolButtonTextUnderIcon)


###################################################
### code chunk number 445: qt-app-statusbar
###################################################
statusbar <- Qt$QStatusBar()
main_window$setStatusBar(statusbar)


###################################################
### code chunk number 446: qt-app-statusbar-temporary
###################################################
statusbar$showMessage("Load complete", 1000)


###################################################
### code chunk number 447: qt-app-statusbar-perm
###################################################
statusbar$addWidget(Qt$QLabel("Ready"))
statusbar$addPermanentWidget(Qt$QLabel("Version 1.0"))


###################################################
### code chunk number 448: qt-app-dockwidget
###################################################
library(qtutils)
device <- QT()
dock <- Qt$QDockWidget("Device 1")
dock$setWidget(device)


###################################################
### code chunk number 449: qt-app-dock-features
###################################################
dock$features <- Qt$QDockWidget$DockWidgetMovable |
                 Qt$QDockWidget$DockWidgetFloatable


###################################################
### code chunk number 450: qt-app-dock-add
###################################################
main_window$addDockWidget(Qt$Qt$LeftDockWidgetArea, dock)


###################################################
### code chunk number 451: qt-app-dock-tabify
###################################################
device2 <- QT()
dock2 <- Qt$QDockWidget("Device 2", device2)
main_window$tabifyDockWidget(dock, dock2)


###################################################
### code chunk number 452: qt-app-floating
###################################################
dock2$floating <- TRUE


