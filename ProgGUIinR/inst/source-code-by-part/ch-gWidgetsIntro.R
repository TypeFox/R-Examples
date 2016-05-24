### R code from vignette source 'ch-gWidgetsIntro.Rnw'

###################################################
### code chunk number 1: ch-gWidgetsIntro.Rnw:17-19
###################################################
options(prompt=" ",continue=" ") 
source("../booktabs.R")


###################################################
### code chunk number 2: Overview.Rnw:58-92 (eval = FALSE)
###################################################
## ## show (linux, mac, windows) x (RGtk2, tcltk, rJava)
## ## also test-33.html in rpad dddirectory for gWidgetrs
## f <- function(os,toolkit,parent = NULL) {
##   w <- gwindow(paste(os,toolkit,sep = ":"), width = 200, parent=parent)
##   g <- ggroup(horiz = FALSE, cont = w)
##   
##   lst <- list(type = "fieldset",
##   label = "argument",
##   children = list(
##   list(type = "gcombobox",
##   label = "combo",
##   items = letters),
##   list(type = "gslider",
##   label = "slider"),
##   list(type = "gedit",
##   label = "edit",
##   text = "edit this")
##   )
##   )
## 
##   gformlayout(lst, cont = g)
##   
##   bg <- ggroup(cont = g)
##   addSpring(bg)
##   gbutton("ok", cont = bg)
## }
## 
## library(gWidgets)
## os <- "Linux"
## for(j in c("RGtk2","tcltk","Qt")) {
##   options(guiToolkit = j)
##   f(os,j)
## }
## 


###################################################
### code chunk number 3: Overview.Rnw:183-217
###################################################
require(gWidgets)
options(guiToolkit="RGtk2")
## 
window <- gwindow("File search", visible=FALSE)
paned <- gpanedgroup(cont = window)
## label and file selection widget
group <- ggroup(cont = paned, horizontal = FALSE)
glabel("Search for (filename):", cont=group, anchor=c(-1,0))
txt_pattern <- gedit("", initial.msg = "Possibly wildcards", 
                    cont = group)
##
glabel("Search in:", cont = group, anchor = c(-1,0))
start_dir <- gfilebrowse(text = "Select a directory ...",
                        quote = FALSE,
                        type = "selectdir", cont = group)
## A button to initiate the search
search_button <- gbutton("Search", cont = group)
addSpring(group)
## Area for output
frame <- gframe("Output:", cont = paned, horizontal = FALSE)
search_results <- gtext("", cont = frame, expand = TRUE)
size(search_results) <- c(350, 200)
## add interactivity
addHandlerChanged(search_button, handler = function(h,...) {
  pattern <- glob2rx(svalue(txt_pattern))
  file_names <- dir(svalue(start_dir), pattern, 
                    recursive = TRUE)
  if(length(file_names))
    svalue(search_results) <- file_names
  else
    galert("No matching files found", parent = window)
})
## display GUI
visible(window) <- TRUE


###################################################
### code chunk number 4: Overview.Rnw:235-237 (eval = FALSE)
###################################################
## gname(some_arguments, handler = NULL, action = NULL, 
##       container = NULL, ...,  toolkit = guiToolkit())


###################################################
### code chunk number 5: Overview.Rnw:296-297 (eval = FALSE)
###################################################
## options(guiToolkit = "RGtk2")


###################################################
### code chunk number 6: Overview.Rnw:386-387 (eval = FALSE)
###################################################
## svalue(txt_pattern)


###################################################
### code chunk number 7: Overview.Rnw:396-397 (eval = FALSE)
###################################################
## svalue(search_results) <- file_names


###################################################
### code chunk number 8: Overview.Rnw:415-416 (eval = FALSE)
###################################################
## visible(window) <- TRUE


###################################################
### code chunk number 9: methods
###################################################
df <- rbind(
            c("\\meth{svalue, svalue\\ASSIGN}", "Get or set widget's main property"),
            c("\\meth{size\\ASSIGN}","Set preferred size request of widget in pixels"),
            c("\\meth{show}", "Show widget if not visible"),
            c("\\meth{dispose}","Destroy widget or its parent"),
            c("\\meth{enabled, enabled\\ASSIGN}","Adjust sensitivity to user input"),
            c("\\meth{visible, visible\\ASSIGN}","Show or hide object or part of object"),
            c("\\meth{focus\\ASSIGN}","Set focus to widget"),
            c("\\meth{insert}","Insert text into a multiline text widget"),
            c("\\meth{font\\ASSIGN}","Set a widget's font"),
            c("\\meth{update}","Update widget value"),
            c("\\meth{isExtant}","Does \\R\\/ object refer to GUI object that still exists"),
            c("",""),
            c("\\meth{[, [\\ASSIGN}","Refers to values in data store"),
            c("\\meth{length}", "\\meth{length} of data store"),
            c("\\meth{dim}","\\meth{dim} of data store"),
            c("\\meth{names}","\\meth{names} of data store "),
            c("\\meth{dimnames}","\\meth{dimnames} of data store"),
            c("",""),
##            c("\\meth{tag, tag\\ASSIGN}","Set an attribute for a widget that persists through copies. Use sparingly."),
##            c("\\meth{defaultWidget\\ASSIGN}","Set widget to have initial focus in a dialog"),
            c("\\meth{getToolkitWidget}","Return underlying toolkit widget for low-level use")
            )
dimnames(df) <- list(NULL,c("Method", "Description"))
cat(booktabs(df, colTypes=c("l","p{0.6\\textwidth}"),
             caption="Generic functions provided or used in the \\pkg{gWidgets} API.",
             label="tab:gWidgets-methods"))


###################################################
### code chunk number 10: Overview.Rnw:531-532 (eval = FALSE)
###################################################
## search_button <- gbutton("Search", cont = group)


###################################################
### code chunk number 11: Overview.Rnw:542-550 (eval = FALSE)
###################################################
## addHandlerChanged(search_button, handler = function(h,...) {
##   pattern <- glob2rx(svalue(txt_pattern))
##   file_names <- dir(svalue(start_dir), pattern,recursive=TRUE)
##   if(length(file_names))
##     svalue(search_results) <- file_names
##   else
##     galert("No matching files found", parent = window)
## })


###################################################
### code chunk number 12: Overview.Rnw:559-560 (eval = FALSE)
###################################################
## prop <- svalue(h$obj)


###################################################
### code chunk number 13: Overview.Rnw:598-611 (eval = FALSE)
###################################################
## search_button <- gbutton("Search", cont = group,
##              handler = function(h,...) {
##                pattern <- glob2rx(svalue(h$action$txt))
##                file_names <- dir(svalue(h$action$dir), 
##                              pattern, recursive = TRUE)
##                if(length(file_names))
##                  svalue(h$action$results) <- file_names
##                else
##                  galert("No matching files found", parent = w)
##              },
##              action = list(txt = txt_pattern, dir = start_dir,
##                results = search_results)
##              )


###################################################
### code chunk number 14: addHandlerXXX
###################################################
df <- rbind(
            c("\\meth{addHandlerChanged}", paste("Primary handler call for when a widget's value is \"changed.\" The interpretation of \"change\" depends on the widget.",sep="")),
            c("\\meth{addHandlerClicked}", paste("Set handler for when widget is clicked with (left)",
                                                 c("mouse button. May return position of click through components"),
                                                 c("\\code{x} and \\code{y} of the \\code{h}-list. "),sep=" ")),
            c("\\meth{addHandlerDoubleclick}","Set handler for when widget is double-clicked."),
            c("\\meth{addHandlerRightclick}","Set handler for when widget is right-clicked."),
            c("\\meth{addHandlerKeystroke}", paste("Set handler for when key is",
                                                   c("pressed. The \\code{key} component is set to this value,"),
                                                   c("if possible."),sep=" ")),
            c("\\meth{addHandlerFocus}","Set handler for when widget gets focus."),
            c("\\meth{addHandlerBlur}","Set handler for when widget loses focus."),
            c("\\meth{addHandlerExpose}","Set handler for when widget is first drawn."),
            c("\\meth{addHandlerUnrealize}","Set handler for when widget is undrawn on screen."),
            c("\\meth{addHandlerDestroy}","Set handler for when widget is destroyed."),
            c("\\meth{addHandlerMouseMotion}","Set handler for when widget has mouse go over it."),
            c("\\meth{addDropSource}","Specify a widget as a drop source."),
            c("\\meth{addDropMotion}","Set handler to be called when an item is dragged over the widget."),
            c("\\meth{addDropTarget}","Set handler to be called on a drop event. Adds the component \\code{dropdata}."),
            c("\\meth{addHandler}","(Not cross-toolkit) Allows one to specify an underlying signal from the graphical toolkit and handler."), 
            c("",""),
            c("\\meth{removeHandler}","Remove a handler from a widget."),
            c("\\meth{blockHandler}","Temporarily block a handler from being called."),
            c("\\meth{unblockHandler}","Restore handler that has been blocked."),
            c("\\meth{addHandlerIdle}","Call a handler during idle time."),
            c("",""),
            c("\\meth{addPopupmenu}","Bind pop-up menu to widget."),
            c("\\meth{add3rdMousePopupmenu}","Bind popup menu to right mouse click.")
            )
dimnames(df) <- list(NULL,c("Method","Description"))
cat(booktabs(df, colTypes=c("l","p{0.6\\textwidth}"),
             caption="Generic functions to add callbacks in \\pkg{gWidgets} API.",
             label="tab:gWidgets-callback-methods"))


###################################################
### code chunk number 15: gfile (eval = FALSE)
###################################################
## if(!is.na(f <- gfile())) source(f)


###################################################
### code chunk number 16: dialogs
###################################################
df <- rbind(
            c("\\constructor{gmessage}", "Dialog to show a message."),
            c("\\constructor{galert}", "Unobtrusive (non-modal) dialog to show a message."),
            c("\\constructor{gconfirm}", "Confirmation dialog."),
            c("\\constructor{ginput}", "Dialog allowing user input."),
            c("\\constructor{gbasicdialog}", "Flexible modal dialog."),
            c("\\constructor{gfile}","File and directory selection dialog.")
            )
colnames(df) <- c("Constructor","Description")
cat(booktabs(df,
             colTypes=c("l","p{0.7\\textwidth}"),
             caption="Table of constructors for basic dialogs in \\pkg{gWidgets}.",
             label="tab:gWidgets-basic-dialogs"))



###################################################
### code chunk number 17: gconfirm (eval = FALSE)
###################################################
## gconfirm("Yes or no? Click one.")


###################################################
### code chunk number 18: gmessage (eval = FALSE)
###################################################
## gmessage("Message goes here", title = "example dialog")


###################################################
### code chunk number 19: Overview.Rnw:774-775 (eval = FALSE)
###################################################
## ret <- gconfirm("Really delete file?", icon = "question")


###################################################
### code chunk number 20: ginput (eval = FALSE)
###################################################
## ret <- ginput("Enter your name", icon = "info")
## if(!is.na(ret)) 
##   message("Hello", ret,"\n")


###################################################
### code chunk number 21: installation
###################################################
df <- rbind(
           Windows = c("Windows", "Installed by \\pkg{RGtk2}",
             "Included with \\pkg{qtbase}",
             "In binary install of R"
             ),
            Linux = c("Linux","Standard","Standard","Standard"),
            Mac = c("OS X", "Download binary .pkg",
              "Vendor supplied",
              "In binary install of R")
            )
dimnames(df) <- list(NULL,c("","Gtk+", "Qt", "Tk"))
cat(booktabs(df, colTypes = c("l", rep("p{0.3\\textwidth}",3)),
             caption = "Installation notes for GUI toolkits.",
             label = "tab:gWidgets-installation"))



###################################################
### code chunk number 22: Containers.Rnw:2-4
###################################################
library(gWidgets)
options("guiToolkit"="RGtk2") ## for examples


###################################################
### code chunk number 23: windowTest1
###################################################
window <- gwindow("Our title", visible = TRUE)


###################################################
### code chunk number 24: windowTest1Container
###################################################
label <- glabel("A child label", container = window)


###################################################
### code chunk number 25: window-visible-false
###################################################
window <- gwindow("Title", visible = FALSE)
## perform layout here ...
visible(window) <- TRUE


###################################################
### code chunk number 26: childWindowTest
###################################################
child_window <- gwindow("A child window", parent = window, 
                        width = 200, height = 200)


###################################################
### code chunk number 27: Containers.Rnw:175-186
###################################################
old_options <- options(error = function() {
  if(msg <- geterrmessage() != "")
    galert(msg, parent = window)
  invisible(msg)
})
#
window <- gwindow( "Popup errors", visible = FALSE,
                  handler = function(h, ...) {
                    ## restore old options when gui is closed
                    options(old_options)       
                  })


###################################################
### code chunk number 28: Containers.Rnw:190-194
###################################################
button <- gbutton("Click for error",  cont = window,
                  handler = function(h, ...) {
                    stop("This is an error")
                  })


###################################################
### code chunk number 29: windowUnrealizeMethod
###################################################
window <- gwindow("Close through the window manager")
id <- addHandlerUnrealize(window, handler = function(h,...) {
  !gconfirm("Really close", parent = h$obj)
})


###################################################
### code chunk number 30: Containers.Rnw:224-241
###################################################
df <- rbind(
            c("\\constructor{gwindow}","Creates a top-level window."),
            c("\\constructor{ggroup}","Creates a box-like container."),
            c("\\constructor{gframe}","Creates a box container with a text label."),
            c("\\constructor{gexpandgroup}",paste("Creates a box container with a label and",
                                                  "a trigger to expand/collapse.",sep=" ")),
            c("\\constructor{glayout}","Creates a grid container."),
            c("\\constructor{gpanedgroup}",
              paste("Creates a container for two child widgets",
                    "with a handle to assign allocation of space.",sep=" ")),
            c("\\constructor{gnotebook}","Creates a tabbed notebook container for holding a collection of child widgets.")
            )
colnames(df) <- c("Constructor","Description")
cat(booktabs(df, 
             colTypes = c("l","p{0.6\\textwidth}"),
             caption="Constructors for container objects.",
             label="tab:gWidgets-container-constructors"))


###################################################
### code chunk number 31: collapseFactor
###################################################
collapseFactor <- function(fac, parent = NULL) {
  out <- character()
  window <- 
    gbasicdialog("Collapse factor levels", parent = parent,
                 handler = function(h,...) {
                   new_fac <- relevel_factor$get_value()
                   out <<- factor(new_fac)
                 })
  group <- ggroup(cont = window)
  relevel_factor <- CollapseFactor$new(fac, cont = group)
  visible(window, set = TRUE)
  out
}


###################################################
### code chunk number 32: Containers.Rnw:299-300 (eval = FALSE)
###################################################
## mtcars$am <- collapseFactor(mtcars$am)


###################################################
### code chunk number 33: Containers.Rnw:316-333
###################################################
## List methods for containers
df <- rbind(c("\\meth{add}",paste("Adds a child object to a parent container.",
                                  "Called when a parent container is specified to the \\args{container}",
                                  "argument of the widget constructor, in which case",
                                  "the \\args{...} arguments are passed to",
                                  "this method.", sep=" ")),
            c("\\meth{delete}", "Removes a child object from a parent container."),
            c("\\meth{dispose}", "Destroys container and children."),
            c("\\meth{enabled\\ASSIGN}", "Sets sensitivity of child components."),
            c("\\meth{visible\\ASSIGN}", "Hides or shows child components.")
            )
colnames(df) <- c("Method","Description")
cat(booktabs(df, 
             colTypes=c("l","p{0.6\\textwidth}"),
             caption="Container methods.", 
             label="tab:gWidgets-container-methods"
             ))


###################################################
### code chunk number 34: cancelOk
###################################################
window <- gwindow("Some buttons", visible = FALSE)
group <- ggroup(horizontal = TRUE, cont = window)
cancel_button <- gbutton("cancel", cont = group)
ok_button <- gbutton("ok", cont = group)
visible(window) <- TRUE


###################################################
### code chunk number 35: Containers.Rnw:392-406
###################################################
## ## Example showing expand argument -- it gets implemented differently in the various
## ## toolkits.
## w <- gwindow("Expand", visible = FALSE)
## g <- ggroup(cont = w, horizontal = FALSE)

## gbutton("Expand = FALSE", expand = FALSE, cont = g)
## gbutton("Expand = TRUE, no fill", expand = TRUE, cont = g)
## gbutton("Expand = TRUE, fill = 'both'", expand = TRUE, fill = "both", cont = g)
## visible(w) <- TRUE

## size(w) <- c(300, 200)
## ## For gWidgetsRGtk2 the buttons all fill in "x", but only the second two are allocated space
## ## For tcltk, the first button does not fill in x (no expand -- no fill)
## for Qt, the first button fils, the second only expands in "y" direction


###################################################
### code chunk number 36: Containers.Rnw:461-471
###################################################
window <- gwindow("Add args")
group <- ggroup(horizontal = FALSE, cont = window)
button1 <- gbutton("expand = FALSE", expand = FALSE, cont = group)
gseparator(cont = group)
button2 <- gbutton("expand = TRUE, fill = 'x'", expand = TRUE, fill = "x", cont = group)
gseparator(cont = group)
button3 <- gbutton("expand = TRUE, fill = 'y'", expand = TRUE, fill = "y", cont = group)
gseparator(cont = group)
button4 <- gbutton("expand = FALSE, anchor = c(1,0)", anchor = c(1,0), cont = group)
update(window)


###################################################
### code chunk number 37: Containers.Rnw:505-513
###################################################
window <- gwindow("Some buttons", visible = FALSE)
group <- ggroup(horizontal = TRUE, spacing = 6, cont = window)
help_button <- gbutton("help", cont = group)
addSpring(group)
cancel_button <- gbutton("cancel", cont = group)
addSpace(group, 12)                      # 6 + 12 + 6 pixels
ok_button <- gbutton("ok", cont = group)
visible(window) <- TRUE


###################################################
### code chunk number 38: Containers.Rnw:589-594
###################################################
window <- gwindow("gframe example")
frame <- gframe("gWidgets Examples:", cont = window)
files <- list.files(system.file("Examples","ch-gWidgets", 
                                package = "ProgGUIinR"))
vars <- gtable(files, cont = frame, expand = TRUE)


###################################################
### code chunk number 39: Containers.Rnw:613-621
###################################################
res <- lm(mpg ~ wt, mtcars)
out <- capture.output(summary(res))
##
window <- gwindow("gexpandgroup example", visible = FALSE)
exp_group <- gexpandgroup("Summary", cont = window)
label <- glabel(out, cont = exp_group)
visible(exp_group) <- TRUE                   # display summary
visible(window) <- TRUE


###################################################
### code chunk number 40: glayoutExample
###################################################
window <- gwindow("glayout example", visible = FALSE)
lyt <- glayout(cont = window, spacing = 5)
right <- c(1,0); left <- c(-1,0)
lyt[1,1, anchor = right] <- "name"
lyt[1,2, anchor = left ] <- gedit("George Washington", 
           cont = lyt)
#
lyt[2,1, anchor = right] <- "rank"
lyt[2,2, anchor = left ] <- gedit("General", cont = lyt)
#
lyt[3,1, anchor = right] <- "serial number"
lyt[3,2, anchor = left ] <- gedit("1", cont = lyt)
visible(window) <- TRUE


###################################################
### code chunk number 41: main_table_prop
###################################################
sapply(lyt[,2], svalue)


###################################################
### code chunk number 42: layout_example
###################################################
examples <- system.file("Examples", "ch-gWidgets", 
                 package = "ProgGUIinR")
files <- list.files(examples)
#
window <- gwindow("gpanedgroup example", visible = FALSE)
paned <- gpanedgroup(cont  =  window)
tbl <- gtable(files, cont = paned)           # left side
txt_widget <- gtext("", cont = paned, expand = TRUE) # right
visible(window) <- TRUE
svalue(paned) <- 0.33                        # after drawing


###################################################
### code chunk number 43: tabbed_notebook
###################################################
window <- gwindow("gnotebook example")
notebook <- gnotebook(cont = window)


###################################################
### code chunk number 44: addPage
###################################################
add_a_page <- function(file_name) {
  f <- system.file(file_name, package = "ProgGUIinR")
  gtext(paste(readLines(f), collapse="\n"), 
        cont = notebook, label = file_name)
}
add_a_page("DESCRIPTION")


###################################################
### code chunk number 45: helpPage
###################################################
lyt <- glayout(cont = notebook, horizontal = FALSE, 
               label = "Help")
lyt[1,1] <- gimage("help", dir = "stock", cont = lyt)
lyt[1,2] <- glabel(paste("To add a page:",
           "Click on a file in the left pane, and its contents",
           "are displayed in a notebook page.", sep = "\n"), 
           cont = lyt)


###################################################
### code chunk number 46: set_notebook_page
###################################################
svalue(notebook) <- length(notebook)


###################################################
### code chunk number 47: rm_notebook_page
###################################################
dispose(notebook)


###################################################
### code chunk number 48: Controls.Rnw:2-5
###################################################
library(gWidgets)
options(guiToolkit = "RGtk2")



###################################################
### code chunk number 49: Controls.Rnw:19-51
###################################################
df <- rbind(
            c("\\constructor{glabel}", "A text label."),
            c("\\constructor{gbutton}", "A button to initiate an action."),
            c("\\constructor{gcheckbox}", "A checkbox."),
            c("\\constructor{gcheckboxgroup}", "A group of checkboxes."),
            c("\\constructor{gradio}", "A radio button group."),
            c("\\constructor{gcombobox}", "A drop-down list of values, possibly editable."),
            c("\\constructor{gtable}", "A table (vector or data frame) of values for selection."),
            c("\\constructor{gslider}", "A slider to select from a sequence value."),
            c("\\constructor{gspinbutton}", "A spinbutton to select from a sequence of values."),
            c("\\constructor{gedit}", "Single line of editable text."),
            c("\\constructor{gtext}", "Multiline text edit area."),
            c("\\constructor{ghtml}", "Display text marked up with HTML."),
            c("\\constructor{gdf}", "Data frame viewer and editor."),
            c("\\constructor{gtree}", "A display for hierarchical data."),
            c("\\constructor{gimage}", "A display for icons and images."),
            c("\\constructor{ggraphics}", "A widget containing a graphics device."),
            c("\\constructor{gsvg}", "A widget to display SVG files."),            
            c("\\constructor{gfilebrowse}", "A widget to select a file or directory."),
            c("\\constructor{gcalendar}", "A widget to select a date."),
            c("\\constructor{gaction}", "A reusable definition of an action."),
            c("\\constructor{gmenubar}", "Add a menu bar to a top-level window."),
            c("\\constructor{gtoolbar}", "Add a toolbar to a top-level window."),
            c("\\constructor{gstatusbar}", "Add a status bar to a top-level window."),
            c("\\constructor{gtooltip}", "Add a tooltip to a widget."),
            c("\\constructor{gseparator}", "A widget to display a horizontal or vertical line.")
            )
colnames(df) <- c("Constructor","Description")
cat(booktabs(df,
             colTypes=c("l","p{0.7\\textwidth}"),
             caption="Table of constructors for control widgets in \\pkg{gWidgets}. Most, but not all, are implemented for each toolkit.",
             label="tab:gWidgets-control-widgets"))


###################################################
### code chunk number 50: Controls.Rnw:117-123
###################################################
window <- gwindow("Make a plot")
group <- ggroup(horizontal = FALSE, cont = window)
glabel("... Fill me in ...", cont = group)
button_group <- ggroup(cont = group)
addSpring(button_group)
parButton <- gbutton("par (mfrow) ...", cont = button_group)


###################################################
### code chunk number 51: Controls.Rnw:128-142
###################################################
addHandlerClicked(parButton, handler = function(h,...) {
  child <- gwindow("Set par values for mfrow", parent = window)
  lyt <- glayout(cont = child)
  lyt[1,1, align = c(-1,0)] <- "mfrow: c(nr,nc)"
  lyt[2,1] <- (nr <- gedit(1, cont = lyt))
  lyt[2,2] <- (nc <- gedit(1, cont = lyt))
  lyt[3,2] <- 
    gbutton("ok", cont = lyt, handler = 
            function(h,...) {
              x <- as.numeric(c(svalue(nr), svalue(nc)))
              par(mfrow = x)
              dispose(child)
            })
})


###################################################
### code chunk number 52: Controls.Rnw:218-228
###################################################
window <- gwindow("label example")
frame <- gframe("Summary statistics:", cont = window)
lyt <- glayout(cont = frame)
lyt[1,1] <- glabel("xbar:", cont = lyt)
lyt[1,2] <- gedit("", cont = lyt)
lyt[2,1] <- glabel("s:", cont = lyt)
lyt[2,2] <- gedit("", cont = lyt)
sapply(lyt[,1], function(i) {
  font(i) <- c(weight = "bold", color = "blue")
})


###################################################
### code chunk number 53: Controls.Rnw:299-306
###################################################
f <- tempfile()
png(f)                                  # not gWidgetstcltk!
hist(rnorm(100))
dev.off()
#
window <- gwindow("Example to show a graphic")
gimage(basename(f), dirname(f), cont = window)


###################################################
### code chunk number 54: ex-gWidgets-add-icons.Rnw:18-21
###################################################
some_colors <- c("black", "red", "blue", "brown",
                "green", "yellow", "purple",
                paste("grey", seq.int(10,90,by=10), sep = ""))


###################################################
### code chunk number 55: ex-gWidgets-add-icons.Rnw:26-37
###################################################
require(grid)
icon_dir <- tempdir(); iconSize <- 16;
make_color_icon <- function(i) {
  filename <- file.path(icon_dir, 
                        sprintf("color-%s.png", i))
  png(file = filename, width = iconSize, height = iconSize)
  grid.newpage()
  grid.draw(rectGrob(gp = gpar(fill = i)))
  dev.off()
  return(filename)
}


###################################################
### code chunk number 56: ex-gWidgets-add-icons.Rnw:43-46
###################################################
icons <- sapply(some_colors, make_color_icon)
icon_names <- sprintf("color-%s", some_colors)
addStockIcons(icon_names, icons)


###################################################
### code chunk number 57: ex-gWidgets-add-icons.Rnw:53-64
###################################################
window <- gwindow("Icon example", visible=FALSE)
callback <- function(h,...) galert(h$action, parent = window)
lyt <- glayout(cont = window, spacing = 0)
for(i in 1:4) {
  for(j in 1:4) {
    ind <- (i - 1) * 4 + j
    lyt[i,j] <- gimage(icons[ind], handler = callback, 
                       action = icon_names[ind], cont = lyt)
  }
}
visible(window) <- TRUE


###################################################
### code chunk number 58: gedit
###################################################
window <- gwindow("Simple gedit example", visible = FALSE)
group <- ggroup(cont = window)
entry <- gedit("", initial.msg = "Enter your name...", 
               cont = group)
visible(window) <- TRUE


###################################################
### code chunk number 59: Controls.Rnw:390-396
###################################################
window <- gwindow("gedit example", visible = FALSE) 
group <- ggroup(cont = window)
glabel("State name:", cont = group)
entry <- gedit("", cont = group)
entry[] <- state.name
visible(window) <- TRUE


###################################################
### code chunk number 60: validationExample
###################################################
require(gWidgets)


###################################################
### code chunk number 61: ex-gWidgets-gedit-validation.Rnw:14-18
###################################################
window <- gwindow("Validation example")
lyt <- glayout(cont = window)
lyt[1,1] <- "R expression:"
lyt[1,2] <- (entry <- gedit("", cont = lyt))


###################################################
### code chunk number 62: ex-gWidgets-gedit-validation.Rnw:33-38
###################################################
require(evaluate)
isValid <- function(e) {
  out <- try(evaluate:::evaluate(e), silent=TRUE)
  !(inherits(out, "try-error") ||  is(out[[2]], "error"))
}


###################################################
### code chunk number 63: validate
###################################################
addHandlerChanged(entry, handler = function(h,...) {
  cur_val <- svalue(entry)
  if(isValid(cur_val)) {
    font(entry) <- c(color = "black")
  } else {
    font(entry) <- c(color = "red")
  }
})


###################################################
### code chunk number 64: Controls.Rnw:445-458
###################################################
df <- rbind(
            c("weight","light, normal, bold"),
            c("style", "normal, oblique, italic"),
            c("family", "normal, sans, serif, monospace"),
            c("size", "a point size, such as 12"),
            c("color", "a named color")
            )
colnames(df) <- c("Attribute","Possible value")
cat(booktabs(df, 
             colTypes=c("l","p{0.6\\textwidth}"),
             caption="Possible specifications for setting font properties. Font values of an object are changed with named lists, as in \\code{font(obj)\\ASSIGN list(weight=\"bold\", size=12, color=\"red\")}.", 
             label="tab:gWidgets-font-properties"
             ))


###################################################
### code chunk number 65: ex-gWidgets-calculator.Rnw:2-4
###################################################
## A calculator layout with gWidgets
library(gWidgets)


###################################################
### code chunk number 66: ex-gWidgets-calculator.Rnw:27-52
###################################################
buttons <- rbind(c(7:9, "(", ")"),
                 c(4:6, "*", "/"),
                 c(1:3, "+", "-"))
#
window <- gwindow("glayout for a calculator", visible = FALSE)
group <- ggroup(cont = window, expand = TRUE, horiz = FALSE)
lyt <- glayout(cont = group, spacing = 2)
                                        
lyt[1, 1:5, anchor = c(-1,0)] <-          # span 5 columns
  (eqn_area <- gedit("", cont = lyt))
lyt[2, 1:5, anchor = c(1,0)] <- 
  (output_area <- glabel("", cont = lyt))
#
button_list <- list()
for(i in 3:5) {
  for(j in 1:5) {
    val <- buttons[i-2, j]
    lyt[i,j] <- (button_list[[val]] <- gbutton(val, cont=lyt))
  }
}
lyt[6,2] <- (button_list[["0"]] <- gbutton("0", cont = lyt))
lyt[6,3] <- (button_list[["."]] <- gbutton(".", cont = lyt))
lyt[6,4:5] <- (eq_button <- gbutton(" = ", cont = lyt))
#
visible(window) <- TRUE


###################################################
### code chunk number 67: ex-gWidgets-calculator.Rnw:60-67
###################################################
add_button <- function(h, ...) {
  cur_expr <- svalue(eqn_area)
  new_char <- svalue(h$obj)              # the button's value
  svalue(eqn_area) <- paste(cur_expr, new_char, sep = "")
  svalue(output_area) <- ""              # clear label 
}
sapply(button_list, addHandlerChanged, handler = add_button)


###################################################
### code chunk number 68: ex-gWidgets-calculator.Rnw:73-88
###################################################
require(evaluate)
addHandlerClicked(eq_button, handler = function(h,...) {
  curExpr <- svalue(eqn_area)
  out <- try(evaluate:::evaluate(curExpr), silent = TRUE)
  if(inherits(out, "try-error")) {
    galert("Parse error", parent = eq_button)
  } else if(is(out[[2]], "error")) {
    msg <- sprintf("Error: %s", out[[2]]$message)
    galert(msg, parent = eq_button)
  } else {
    svalue(output_area) <- out[[2]]
    svalue(eqn_area) <- ""            # restart
  }
})
                  


###################################################
### code chunk number 69: Controls.Rnw:532-535
###################################################
window <- gwindow("Checkbox example with toggle button")
check_box <- gcheckbox("Thresh", checked = TRUE, 
                       use.togglebutton = TRUE, cont = window)


###################################################
### code chunk number 70: Controls.Rnw:549-555
###################################################
window <- gwindow("checkbox example")
check_button <- gcheckbox("label", cont = window, 
                          handler = function(h,...) {
                            if(svalue(h$obj)) # it is checked
                              print("define handler here")
                          })


###################################################
### code chunk number 71: Controls.Rnw:575-578
###################################################
window <- gwindow("Radio button example")
radio_button <- gradio(c("Color", "Grayscale"), selected = 2, 
                       horizontal = FALSE, cont = window)


###################################################
### code chunk number 72: Controls.Rnw:623-628
###################################################
window <- gwindow("Checkbox group example")
check_box_group <-
  gcheckboxgroup(c("Flip","Flop"), horizontal = FALSE, 
                 checked = c(FALSE, TRUE), cont = window)



###################################################
### code chunk number 73: Controls.Rnw:647-650
###################################################
svalue(check_box_group) <- c("Flop")
svalue(check_box_group) <- c(FALSE, TRUE)
svalue(check_box_group, index = TRUE) <- 2


###################################################
### code chunk number 74: Controls.Rnw:710-712
###################################################
window <- gwindow("gcombobox example")
combo_box <- gcombobox(c("None","Low","High"), cont = window)


###################################################
### code chunk number 75: comboboxExample
###################################################
nms <- getStockIcons()                  # gWidgets icons
DF <- data.frame(names = names(nms), icons = names(nms), 
                 stringsAsFactors = FALSE)
window <- gwindow("Combo box with icons example")
combo_box <- gcombobox(DF, cont = window)


###################################################
### code chunk number 76: Controls.Rnw:772-775
###################################################
avail_DFs <- function() {
  c("", ".GlobalEnv", ProgGUIinR:::avail_dfs(.GlobalEnv))
}


###################################################
### code chunk number 77: Controls.Rnw:777-781
###################################################
##' get numeric variables
##'
##' @param where ".GlobalEnv" or name of data frame in global workspace
##' @return vector of numeric variable names


###################################################
### code chunk number 78: getNumeric
###################################################
get_numeric <- function(where) {
  val <- get(where, envir = .GlobalEnv)
  ProgGUIinR:::find_vars(val, is.numeric)
}


###################################################
### code chunk number 79: Controls.Rnw:791-817
###################################################
window <- gwindow("Find the mean", visible = FALSE)
group <- ggroup(cont = window, horizontal = FALSE)
group1 <- ggroup(cont = group)
glabel("Select data frame:", cont = group1)
df_combo_box <- gcombobox(avail_DFs(), cont = group1)
##
frame <- gframe("Arguments:", cont = group, horizontal=FALSE)
enabled(frame) <- FALSE
lyt <- glayout(cont = frame, expand = TRUE)
widget_list <- list() 
##
lyt[1,1] <- "x"
lyt[1,2] <- (widget_list$x <- gcombobox("           ",
                                        cont = lyt))
##
lyt[2,1] <- "trim"
lyt[2,2] <- 
  (widget_list$trim <- gslider(from = 0, to = 0.5, by = 0.01,
                               cont = lyt))
##
lyt[3,1] <- "na.rm"
lyt[3,2] <- 
  (widget_list$na.rm <- gcheckbox("", checked = TRUE, 
                                  cont = lyt))
group2 <- ggroup(cont = group)
compute_button <- gbutton("compute", cont = group2)


###################################################
### code chunk number 80: Controls.Rnw:829-837
###################################################
addHandlerChanged(df_combo_box, handler = function(h,...) {
  val <- svalue(h$obj)
  enabled(frame) <- val !=""
  enabled(compute_button) <- val != ""
  if(val != "") 
    widget_list$x[] <- get_numeric(val)
  svalue(widget_list$x, index = TRUE) <- 0
})


###################################################
### code chunk number 81: computeHandler
###################################################
addHandlerChanged(compute_button, handler = function(h,...) {
  out <- lapply(widget_list, svalue)
  out$x <- get(out$x, get(svalue(df_combo_box),
                          envir = .GlobalEnv))
  print(do.call(mean.default, out))
})


###################################################
### code chunk number 82: visible
###################################################
visible(window) <- TRUE


###################################################
### code chunk number 83: brightness
###################################################
window <- gwindow("Slider example")
brightness <- gslider(from = -1, to = 1, by = .05, value = 0, 
   handler = function(h,...) {
     cat("Update picture with brightness", svalue(h$obj),"\n")
   }, cont = window)


###################################################
### code chunk number 84: slider
###################################################
window <- gwindow("Add a label to the slider", visible=FALSE)
group <- ggroup(cont = window, expand = TRUE)
slider <- gslider(from = 0, to = 100, by = 1, cont = group, 
                  expand = TRUE)
label <- glabel(sprintf("%3d", svalue(slider)), cont = group)
font(label) <- c(family = "monospace")
addHandlerChanged(slider, function(h,...) {
  svalue(h$action) <- sprintf("%3d", svalue(h$obj))
  }, action = label)
visible(window) <- TRUE


###################################################
### code chunk number 85: spinbutton
###################################################
window <- gwindow("Spin button example")
spin_button <- gspinbutton(from = 0, to = 10, by = .05, 
                           value = 1, cont = window)


###################################################
### code chunk number 86: readImage
###################################################
## stub for readImage from EBImage
readImage <- function() {}


###################################################
### code chunk number 87: openImage (eval = FALSE)
###################################################
## f <- gfile("Open an image file",
##            type="open",
##            filter=list("Image file" = list(
##                        patterns = c("*.gif", "*.jpg", "*.png")
##                        ),
##              "All files" = list(patterns = c("*"))
##              ))
## if(!is.na(f)) 
##   readImage(f) ## ...


###################################################
### code chunk number 88: advSearch
###################################################
## thanks to Richie Cotton for this example
window <- gwindow("File search", visible = FALSE)

paned <- gpanedgroup(cont = window)

nested_group <- ggroup(cont = paned, horizontal = FALSE)
glabel("Search for (filename):", cont = nested_group, anchor = c(-1,0))
txt_pattern <- gedit("", initial.msg = "Possibly wildcards", cont = nested_group)

glabel("Search in:", cont = nested_group, anchor = c(-1,0))
start_dir <- gfilebrowse(text = "Select a directory ...",
                        quote = FALSE,
                        type = "selectdir", cont = nested_group)


###################################################
### code chunk number 89: advSearch
###################################################
adv_search <- gexpandgroup("Advanced search:", 
                           cont = nested_group)
visible(adv_search) <- FALSE
lyt <- glayout(cont = adv_search)
lyt[1,1] <- "Recursive"
lyt[1,2] <- (adv_rec <- 
   gcheckbox("search directories", checked = TRUE, cont=lyt))
lyt[2,1] <- "Size"
lyt[2,2] <- (adv_size <- 
   gcombobox(c("", "small", "medium", "large"),  cont = lyt))
lyt[3,1] <- "All files"
lyt[3,2] <- (adv_visible <- 
   gradio(c(TRUE, FALSE), horizontal = TRUE, cont = lyt))
lyt[4,1] <- "Last modified"
lyt[4,2] <- (adv_modified <- 
             gcalendar("", format = "%Y-%m-%d", cont = lyt))


###################################################
### code chunk number 90: ex-gWidgets-file-search.Rnw:62-67
###################################################
search_btn <- gbutton("Search", cont = nested_group)
addSpring(nested_group)

frame <- gframe("Output:", cont = paned, horizontal = FALSE)
search_results <- gtext("", cont = frame, expand = TRUE, fill = "both")


###################################################
### code chunk number 91: ex-gWidgets-file-search.Rnw:74-96
###################################################
addHandlerChanged(search_btn, handler=function(h,...) {
  pattern <- glob2rx(svalue(txt_pattern))
  start_at <- svalue(start_dir)
  modified <- NULL
  size <- NULL

  ## from advanced
  subfolders <- svalue(adv_rec)
  visible <- svalue(adv_visible)
  if((tmp <- svalue(adv_size)) != "") size <- tmp
  if(!is.na(tmp <- svalue(adv_modified))) modified <- tmp
  
  ## function call
  file_names <- file_search(pattern, start_at, subfolders, 
                        modified = modified,
                        size = size, visible = visible)
  dispose(search_results)                # clear
  if(length(file_names))
    svalue(search_results) <- file_names
  else
    galert("No matching files found", parent = window)
})


###################################################
### code chunk number 92: ex-gWidgets-file-search.Rnw:99-100
###################################################
visible(window) <- TRUE


###################################################
### code chunk number 93: file_search
###################################################
#If you want to return to this example later, here's a function of mine that implements more advanced searching.  (This basically mimics the Windows XP search functionality.)
#I thought this might be good for showing off gexpandgroup (to show/hide the advanced options) and gcalendar.
file_search <- function(pattern, start_dir = getwd(), subfolders = TRUE,
                        modified,
                        modified_format  =  "%Y-%m-%d",
                        size,
                        visible = TRUE)
{
  file_names <- dir(
    start_dir, 
    pattern, 
    recursive = subfolders, 
    all.files = visible
  )
  filter_by_date <- !missing(modified) && !is.null(modified)
  filter_by_size <- !missing(size) && !is.null(size)
  if(filter_by_date || filter_by_size)
  {
    finfo <- file.info(file_names)     
  }
  if(filter_by_date)
  {
    now <- Sys.time()
    one_day <- 60 * 60 * 24
    time_range <- switch(modified,
      week  = c(now - 7 * one_day, now),
      month = c(now - 31 * one_day, now),
      year  = c(now - 365 * one_day, now),
      {
        modified <- strptime(modified, modified_format)
        if(length(modified) == 1L) c(modified, now) else modified
      }
    )
    file_names <- file_names[finfo$mtime > time_range[1] && finfo$mtime < time_range[2]]
  }
  if(filter_by_size)
  {
    size <- switch(size,
      small  = c(0, 100 * 1024),
      medium = c(100 * 1024, 1024 * 1024),
      large  = c(1024 * 1024, Inf),
      size
    )
    file_names <- file_names[finfo$size > size[1] && finfo$size < size[2]]  
  }
  file_names
}



###################################################
### code chunk number 94: Controls.Rnw:1117-1119
###################################################
window <- gwindow("gtable example")
DFs <- gtable(ProgGUIinR:::avail_dfs(), cont = window)


###################################################
### code chunk number 95: updateDfs
###################################################
updateDfs <- function() {
  DFs[] <- ProgGUIinR:::avail_dfs()
}


###################################################
### code chunk number 96: Controls.Rnw:1189-1194
###################################################
addHandlerDoubleclick(DFs, handler = function(h,...) {
  val <- svalue(h$obj)
  print(summary(get(val, envir = .GlobalEnv))) # some action

})


###################################################
### code chunk number 97: initialize
###################################################
initialize <- function(fac, cont = gwindow()) {
  old <<- as.character(fac)
  make_gui(cont)
  callSuper()
}


###################################################
### code chunk number 98: ex-gWidgets-collapse-factor.Rnw:41-78
###################################################
make_gui <- function(cont) {
  group <- gpanedgroup(cont = cont)
  levs <- sort(unique(as.character(old)))
  DF <- data.frame(original = levs,
                   new = levs, stringsAsFactors = FALSE)
  #
  widget <<- tbl <- gtable(DF, cont = group,  multiple = TRUE)
  size(tbl) <- c(300, 200)
  #
  nested_group <- ggroup(cont = group, horizontal = FALSE)
  instructions <- gettext("Select levels, then\ 
enter a new combined level
by typing or selecting a level and then enter")
  #
  glabel(instructions, cont = nested_group)
  combo_box <- gcombobox(levs, selected = 0, editable = TRUE, 
                         cont = nested_group)
  enabled(combo_box) <- FALSE
  #
  addHandlerClicked(widget, function(h,...) {
    ind <- svalue(widget, index = TRUE)
    enabled(combo_box) <- (length(ind) > 0)
  })
  ##
  addHandlerChanged(combo_box, handler = function(h,...) {
    ind <- svalue(tbl, index = TRUE)
    if(length(ind) == 0) 
      return()
    #
    tbl[ind,2] <- svalue(combo_box)
    svalue(tbl, index = TRUE) <- 0
    blockHandler(combo_box)
    combo_box[] <- sort(unique(tbl[,2]))
    svalue(combo_box, index = TRUE) <- 0
    unblockHandler(combo_box)
  })
}


###################################################
### code chunk number 99: ex-gWidgets-collapse-factor.Rnw:83-92
###################################################
get_value <- function() {
  "Return factor with new levels"
  old_levels <- widget[,1]
  new_levels <- widget[,2]
  new <- old
  for(i in seq_along(old_levels)) # one pass
    new[new == old_levels[i]] <- new_levels[i]
  factor(new)
}


###################################################
### code chunk number 100: ex-gWidgets-collapse-factor.Rnw:97-107
###################################################
CollapseFactor <- setRefClass("CollapseFactor",
                              fields = list(
                               old = "ANY",
                               widget = "ANY"
                               ),
                             methods = list(
                               initialize = initialize,
                               make_gui = make_gui,
                               get_value = get_value
                             ))


###################################################
### code chunk number 101: ex-gWidgets-collapse-factor.Rnw:111-132 (eval = FALSE)
###################################################
## ##' Collapse a factor
## ##'
## ##' Collapse a factor through a GUI
## ##' @param f factor to collapse
## ##' @param parent optional, where to place editing dialog
## ##' @return releved factor
## collapseFactor <- function(f, parent = NULL) {
##   out <- character()
##   w <- gbasicdialog("Collapse factor levels", parent = parent,
##                     handler = function(h,...) {
##                       new_f <- relf$get_value()
##                       out <<- factor(new_f)
##                     })
##   g <- ggroup(cont = w)
##   relf <- CollapseFactor$new(f, cont = g)
##   visible(w, set = TRUE)
##   out
## }
##   
## ## test it out
## mtcars$am <- collapseFactor(mtcars$am)


###################################################
### code chunk number 102: Controls.Rnw:1236-1240
###################################################
require(MASS)
window <- gwindow("gtable example")
tbl <- gtable(Cars93, chosencol = 1, filter.column = 3, 
              cont = window)


###################################################
### code chunk number 103: Controls.Rnw:1246-1250
###################################################
addHandlerChanged(tbl, handler = function(h,...) {
  val <- svalue(h$obj, drop = FALSE)
  cat(sprintf("You selected the %s %s", val[,1], val[,2]))
})


###################################################
### code chunk number 104: ex-gWidgets-filter-gtable.Rnw:17-18
###################################################
library(gWidgets)


###################################################
### code chunk number 105: ex-gWidgets-filter-gtable.Rnw:24-25
###################################################
options(repos="http://streaming.stat.iastate.edu/CRAN")


###################################################
### code chunk number 106: ex-gWidgets-filter-gtable.Rnw:27-28
###################################################
avail_pkgs <- available.packages()       # pick a cran site


###################################################
### code chunk number 107: ex-gWidgets-filter-gtable.Rnw:32-37
###################################################
window <- gwindow("test of filter")
group <- ggroup(cont = window, horizontal = FALSE)
entry <- gedit("", cont = group)
tbl <- gtable(avail_pkgs, cont = group, filter.FUN = "manual",
              expand = TRUE)


###################################################
### code chunk number 108: ex-gWidgets-filter-gtable.Rnw:49-52
###################################################
our_match <- function(cur_val, vals) {
  grepl(cur_val, vals)
}


###################################################
### code chunk number 109: ex-gWidgets-filter-gtable.Rnw:61-66
###################################################
id <- addHandlerKeystroke(entry, handler = function(h, ...) {
  vals <- tbl[, 1, drop = TRUE]
  cur_val <- svalue(h$obj)
  visible(tbl) <- our_match(cur_val, vals)
})


###################################################
### code chunk number 110: ex-gWidgets-ws-browser.Rnw:50-63 (eval = FALSE)
###################################################
## ## Code tomake figure
## require(iv)
## observer = umlclass("Observer", 
##   c("update(...)"),
##   c("o"))
## observable = umlclass("Observable", 
##   c("add_observer(o)", "remove_observer(o)", "notify_observers(...)"), 
##   "observers")
## m <- uml(observable, observer)
## umlgui(m)
## ##
## umlpng(m)
## system("mv model.png fig-gWidgets-observer-observable-uml.png")


###################################################
### code chunk number 111: WSModel
###################################################
library(objectProperties)
require(digest)
WSModel <- setRefClass("WSModel",
              fields = c(
                properties(list(ws_objects = "character")),
                ws_objects_digests = "character"
                ))


###################################################
### code chunk number 112: initialize_refresh
###################################################
WSModel$methods(
       .get_objects_digests = function() {
         "Helper function to return list with names, digests"
         items <- ls(envir = .GlobalEnv)
         objects <- mget(items, .GlobalEnv)
         trim <- !sapply(objects, is, class2 = "refClass")
         list(items[trim],
              sapply(objects[trim], digest))
       },
       initialize = function() {
         objs <- .get_objects_digests() # call helper
         initFields(ws_objects = objs[[1]],
                    ws_objects_digests = objs[[2]])
         callSuper()
       },
       refresh = function() {
         objs <- .get_objects_digests()                           
         cur_objects <- objs[[1]]
         cur_digests <- objs[[2]]
         ## changes?
         if(length(cur_digests) != ws_objects_digests ||
            length(ws_objects_digests) == 0 ||
            any(cur_digests != ws_objects_digests)) {
           ws_objects <<- cur_objects # signal
           ws_objects_digests <<- cur_digests
         }})


###################################################
### code chunk number 113: "get_objects"
###################################################
WSModel$methods(
        get = function(klass) {
          "klass a string, such as 'numeric' or '!function'"
          if(missing(klass) || length(klass) == 0)
            return(ws_objects)
          ## if we have klass, more work
          ind <- sapply(mget(ws_objects, .GlobalEnv), 
                        function(x) {
                          any(sapply(klass, function(j)  {
                            if(grepl("^!", j))
                              !is(x, substr(j, 2, nchar(j)))
                            else
                              is(x, j)
                          }))
                        })
          ##
          if(length(ind))
            ws_objects[ind]
          else
            character(0)
        })


###################################################
### code chunk number 114: ex-gWidgets-ws-browser.Rnw:195-199
###################################################
WSModel$methods(
                add_observer = function(FUN, ...) {
                  .self$ws_objectsChanged$connect(FUN, ...)
                })


###################################################
### code chunk number 115: View
###################################################
WSView <- setRefClass("WSView",
                      methods = list(
                        update = function(model) {
                          "Subclass this"
                        },
                        set_model = function(model) {
                          FUN <- function() .self$update(model)
                          model$add_observer(FUN)
                        }
                        ))


###################################################
### code chunk number 116: WidgetView
###################################################
WidgetView <- 
  setRefClass("WidgetView",
              contains = "WSView",
              fields = list(
                klass = "character", # which classes to show
                widget = "ANY"
                ),
              methods = list(
                initialize = function(parent, model, 
                    klass=character(0), ...) {
                  if(!missing(model)) set_model(model)
                  if(!missing(parent))init_widget(parent, ...)
                  initFields(klass=klass)
                  update(model)
                  callSuper()
                },
                init_widget = function(parent, ...) {
                  "Initialize widget"
                }))


###################################################
### code chunk number 117: ex-gWidgets-ws-browser.Rnw:249-251
###################################################
library(gWidgets)
options(guiToolkit="RGtk2")


###################################################
### code chunk number 118: TableView
###################################################
TableView <-
  setRefClass("TableView",
        contains = "WidgetView",
        methods = list(
          init_widget = function(parent, ...) {
            widget <<- gtable(makeDataFrame(character(0)),
                              cont = parent, ...)
          },
          update = function(model, ...) {
            widget[] <<- makeDataFrame(model$get(klass))
          }))


###################################################
### code chunk number 119: size_of
###################################################
size_of <- function(x, ...) UseMethod("size_of")
size_of.default <- function(x, ...) "NA"
size_of.character <- size_of.numeric <- 
  function(x, ...) sprintf("%s elements", length(x))
size_of.matrix <- function(x, ...) 
  sprintf("%s x %s", nrow(x), ncol(x))


###################################################
### code chunk number 120: ex-gWidgets-ws-browser.Rnw:290-298
###################################################
## nicer way to do this.
nsprintf <- function(msg1, msg2, x) {
  sprintf(ngettext(x, msg1, msg2), x)
}
## hidden, too long
size_of.factor <- function(x, ...) nsprintf("Factor with %s level", "Factor with %s levels", length(levels(x)))
size_of.data.frame <- function(x, ...) sprintf("%s x %s", nrow(x), ncol(x))
size_of.list <- function(x, ...) nsprintf("%s component", "%s components", length(x))


###################################################
### code chunk number 121: short_description
###################################################
short_description <- function(x, ...) 
  UseMethod("short_description")
short_description.default <- function(x, ...) "R object"
short_description.numeric <- function(x, ...) "Numeric vector"
short_description.integer <- function(x, ...) "Integer"


###################################################
### code chunk number 122: makeDataFrame
###################################################
makeDataFrame <- function(x, envir = .GlobalEnv) {
  DF <- data.frame(variable = character(0),
                   size = character(0), 
                   description = character(0), 
                   class = character(0),
                   stringsAsFactors = FALSE)
  if(length(x)) {
    l <- mget(x, envir)
    short_class <- function(x) class(x)[1]
    DF <- data.frame(variable = x,
                    size = sapply(l, size_of),
                    description=sapply(l, short_description),
                    class = sapply(l, short_class),
                    stringsAsFactors = FALSE)
  }
  DF
}


###################################################
### code chunk number 123: DfView
###################################################
DfView <-
  setRefClass("DfView",
        contains = "WidgetView",
        methods = list(
          initFields = function(...) klass <<- "data.frame", 
          init_widget = function(parent, ...) {
            DF <- data.frame("Data frames" = character(0),
                            stringsAsFactors = FALSE)
            widget <<- gcombobox(DF, cont = parent, ...)
          },
          update = function(model, ...) {
            widget[] <<- model$get(klass)
          }
          ))


###################################################
### code chunk number 124: testit
###################################################
window <- gwindow()
notebook <- gnotebook(cont = window)
##
model <- WSModel$new()
## basic view of certain classes
view <- TableView$new(parent = notebook, model = model, 
                      label = "data",
                      klass=c("factor","numeric", "character", 
                        "data.frame", "matrix", "list"))
## view of non functions
view1 <- TableView$new(parent = notebook, model = model, 
                       label = "not a function",
                       klass = "!function"
                       )
## view of all
view2 <- TableView$new(parent = notebook, model = model,
                       label = "all")
## a bit contrived here, but useful elsewhere
view3 <- DfView$new(parent = notebook, model = model,
                    label = "data frames")
#
model$refresh()                          
svalue(notebook) <- 1


###################################################
### code chunk number 125: Controls.Rnw:1296-1307
###################################################
offspring <- function(path = character(0), lst, ...) {
  if(length(path))
    obj <- lst[[path]]
  else
      obj <- lst
  #
  f <- function(i) is.recursive(i) && !is.null(names(i))
  data.frame(comps = names(obj), 
             hasOffspring = sapply(obj, f),
             stringsAsFactors = FALSE)
}


###################################################
### code chunk number 126: Controls.Rnw:1317-1321
###################################################
lst <- list(a = "1", b =  list(a = "2", b = "3", 
                       c = list(a = "4")))
window <- gwindow("Tree test")
tree <- gtree(offspring, offspring.data = lst, cont = window)


###################################################
### code chunk number 127: ex-gWidgets-gtree.Rnw:1-2
###################################################
require(gWidgets)


###################################################
### code chunk number 128: party (eval = FALSE)
###################################################
## require(party)
## data("GlaucomaM", package = "ipred")      # load data
## gt <- ctree(Class ~ ., data = GlaucomaM)  # fit model


###################################################
### code chunk number 129: offspring
###################################################
offspring <- function(key, offspring.data) {
  if(missing(key) || length(key) == 0)  # which party node?
    node <- 1
  else
    node <- as.numeric(key[length(key)]) # key is a vector

  if(nodes(gt, node)[[1]]$terminal)    # return if terminal
    return(data.frame(node = node, hasOffspring = FALSE,
                      description = "terminal",
                      stringsAsFactors = FALSE))

  DF <- data.frame(node = integer(2), hasOffspring=logical(2),
                   description = character(2), 
                   stringsAsFactors = FALSE)
  ## party internals
  children <-  c("left","right")
  ineq <- c(" <= "," >  ")
  varName <- nodes(gt, node)[[1]]$psplit$variableName
  splitPoint <- nodes(gt, node)[[1]]$psplit$splitpoint

  for(i in 1:2) {
    DF[i,1] <- nodes(gt, node)[[1]][[children[i]]][[1]]
    DF[i,2] <- !nodes(gt, DF[i,1])[[1]]$terminal
    DF[i,3] <- paste(varName, splitPoint, sep = ineq[i])
  }
  DF                                    # returns a data frame
}


###################################################
### code chunk number 130: makeGUI (eval = FALSE)
###################################################
## window <- gwindow("Example of gtree")
## group <- ggroup(cont = window, horizontal = FALSE)
## label <- glabel("Click on the tree to investigate the
##  partition", cont = group)
## tree <- gtree(offspring, cont = group, expand = TRUE)


###################################################
### code chunk number 131: eventHandler
###################################################
addHandlerDoubleclick(tree, handler = function(h,...) {
  node <- as.numeric(svalue(h$obj))
  if(nodes(gt, node)[[1]]$terminal) {   # if terminal plot
    weights <- as.logical(nodes(gt,node)[[1]]$weights)
    plot(response(gt)[weights, ])
  }})


###################################################
### code chunk number 132: Controls.Rnw:1391-1401
###################################################
window <- gwindow("gaction example")
action <- gaction("click me", tooltip = "Click for a message", 
                  icon = "ok", 
                  handler = function(h, ...) {
                    print("Hello")
                  },
                  parent = window)
button <- gbutton(action = action, cont = window)
## .. to change
enabled(action) <- FALSE                     # can't click now


###################################################
### code chunk number 133: Controls.Rnw:1424-1441
###################################################
stub <- function(h,...) gmessage("called handler", 
                                 parent = window)
action_list = list(
  new = gaction(label = "new", icon = "new", 
    handler = stub, parent = window),
  open = gaction(label = "open", icon = "open", 
    handler = stub, parent = window),
  save = gaction(label = "save", icon = "save", 
    handler = stub, parent = window),
  save.as = gaction(label = "save as...", icon = "save as...", 
    handler = stub, parent = window),
  quit = gaction(label = "quit", icon = "quit", 
    handler = function(...) dispose(window), parent = window),
  cut = gaction(label = "cut", icon = "cut", 
    handler = stub, parent = window)
  )



###################################################
### code chunk number 134: Controls.Rnw:1445-1451
###################################################
window <- gwindow("gtoolbar example")
tool_bar_list<- c(action_list[c("new","save")], 
                 sep = gseparator(), 
                 action_list["quit"])
tool_bar <- gtoolbar(tool_bar_list, cont = window)
gtext("Lorem ipsum ...", cont = window)


###################################################
### code chunk number 135: Controls.Rnw:1478-1491
###################################################
menu_bar_list <- list(file = list(
             new = action_list$new,
             open = action_list$open,
             save = action_list$save,
             "save as..." = action_list$save.as,
             sep = gseparator(),
             quit = action_list$quit
             ),
           edit = list(
             cut = action_list$cut
             )
           )



###################################################
### code chunk number 136: Controls.Rnw:1494-1498
###################################################
window <- gwindow("Menu bar example")
menu_bar <- gmenu(menu_bar_list, cont = window)
tool_bar <- gtoolbar(tool_bar_list, cont = window)
txt_widget <- gtext("", cont = window, expand = TRUE)


###################################################
### code chunk number 137: Controls.Rnw:1527-1535
###################################################
no_changes <- c("save","save.as","cut")
keyhandler <- function(...) {
  for(i in no_changes)
    enabled(action_list[[i]]) <- 
      (nchar(svalue(txt_widget)) > 0)
}
addHandlerKeystroke(txt_widget, handler = keyhandler)
keyhandler()


###################################################
### code chunk number 138: Controls.Rnw:1554-1565
###################################################
window <- gwindow("Popup example")
button <- gbutton("click me or right click me", cont = window, 
                  handler = function(h, ...) {
                    cat("You clicked me\n")
                  })
f <- function(h,...) cat("you right clicked on", h$action)
menu_bar_list <- 
  list(one = gaction("one", action = "one", handler = f),
       two = gaction("two", action = "two", handler = f)
       )
add3rdMousePopupmenu(button, menu_bar_list)


###################################################
### code chunk number 139: CompoundWidgets.Rnw:6-23
###################################################
df <- rbind(
            c("\\constructor{ggraphics}", "Embeddable graphics device"),
            c("\\constructor{ggraphicsnotebook}", "Notebook for multiple devices"),
            c("\\constructor{gdf}", "Data frame editor"),
            c("\\constructor{gdfnotebook}", "Notebook for multiple \\code{gdf} instances"),
            c("\\constructor{gvarbrowser}", "GUI for browsing variables in the workspace"),
            c("\\constructor{gcommandline}", "Command line widget"),
            c("\\constructor{gformlayout}", "Creates a GUI from a list specifying layout"),
            c("\\constructor{ggenericwidget}", "Creates a GUI for a function based on its formal arguments or a defining list")
            )
##             ,c("\\constructor{ghelp}", "GUI for a help page"),
##             c("\\constructor{ghelpbrowser}", "A help browser")
colnames(df) <- c("Constructor","Description")
cat(booktabs(df,
             colTypes=c("l","p{0.7\\textwidth}"),
             caption="Table of constructors for \\R-specific widgets in \\pkg{gWidgets}",
             label="tab:gWidgets-compound-widgets"))


###################################################
### code chunk number 140: ggraphicsExample (eval = FALSE)
###################################################
## library(gWidgets); options(guiToolkit = "RGtk2")
## window <- gwindow("ggraphics example", visible = FALSE)
## plot_device <- ggraphics(cont = window)
## x <- mtcars$wt; y <- mtcars$mpg
## #
## addHandlerClicked(plot_device, handler = function(h, ...) {
##   cat(sprintf("You clicked %.2f x %.2f\n", h$x, h$y))
## })
## addHandlerChanged(plot_device, handler = function(h,...) {
##   rx <- h$x; ry <- h$y
##   if(diff(rx) > diff(range(x))/100 && 
##      diff(ry) > diff(range(y))/100) {
##     ind <- rx[1] <= x & x <= rx[2] & ry[1] <=y & y <= ry[2]
##     if(any(ind))
##       print(cbind(x = x[ind], y = y[ind]))
##   }
## })
## visible(window) <- TRUE
## #
## plot(x, y)


###################################################
### code chunk number 141: tkrplot
###################################################
options(guiToolkit = "tcltk"); require(tkrplot)
window <- gwindow("How to embed tkrplot", visible = FALSE)
group <- ggroup(cont = window, horizontal = FALSE)
bb <- 1
img <- tkrplot(getToolkitWidget(group), 
               fun = function() plot(1:20,(1:20)^bb))
add(group, img)
f <- function(...) {
  b <- svalue(slider)
  print(b)
  if (b != bb) {
    bb <<- b
    tkrreplot(img)
  }
}
slider <- gslider(from = 0.05, to = 2, by = 0.05, cont = group, 
                  handler = f, expand = TRUE)
visible(window) <- TRUE



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


###################################################
### code chunk number 156: gdfExample
###################################################
window <- gwindow("gdf example")
DF <- gdf(mtcars, cont = window)
## ... make some edits ...
new_data_frame <- DF[,]                   # store changes


###################################################
### code chunk number 157: ex-gWidgets-gvarbrowser-dnd.Rnw:1-2
###################################################
library(gWidgets)


###################################################
### code chunk number 158: ex-gWidgets-gvarbrowser-dnd.Rnw:10-19
###################################################
window <- gwindow("Drag-and-drop example")
group <- ggroup(cont = window)
workspace_browser <- gvarbrowser(cont = group)
nested_group <- ggroup(horizontal = FALSE, cont = group, 
                       expand = TRUE)
ggraphics(cont = nested_group)
xlabel <- glabel("", cont = nested_group)
ylabel <- glabel("", cont = nested_group)
clear <- gbutton("clear", cont = nested_group)


###################################################
### code chunk number 159: ex-gWidgets-gvarbrowser-dnd.Rnw:24-31
###################################################
init_txt <- "<Drop %s variable here>"
initUI <- function(...) {
  svalue(xlabel) <- sprintf(init_txt, "x")
  svalue(ylabel) <- sprintf(init_txt, "y")
  enabled(ylabel) <- FALSE
}
initUI()                                # initial call


###################################################
### code chunk number 160: ex-gWidgets-gvarbrowser-dnd.Rnw:35-36
###################################################
addHandlerClicked(clear, handler = initUI)


###################################################
### code chunk number 161: ex-gWidgets-gvarbrowser-dnd.Rnw:44-60
###################################################
updateUI <- function(...) {
  if(grepl(svalue(xlabel), sprintf(init_txt, "x"))) {
    ## none set
    enabled(ylabel) <- FALSE
  } else if(grepl(svalue(ylabel), sprintf(init_txt, "y"))) {
    ## x, not y
    enabled(ylabel) <- TRUE
    x <- eval(parse(text = svalue(xlabel)), envir=.GlobalEnv)
    plot(x, xlab = svalue(xlabel))
  } else {
    enabled(ylabel) <- TRUE    
    x <- eval(parse(text = svalue(xlabel)), envir=.GlobalEnv)
    y <- eval(parse(text = svalue(ylabel)), envir=.GlobalEnv)
    plot(x, y, xlab = svalue(xlabel), ylab = svalue(ylabel))
  }
}


###################################################
### code chunk number 162: ex-gWidgets-gvarbrowser-dnd.Rnw:84-90
###################################################
dropHandler <- function(h,...) {
  svalue(h$obj) <- h$dropdata
  updateUI()
}
addDropTarget(xlabel, handler = dropHandler)
addDropTarget(ylabel, handler = dropHandler)


###################################################
### code chunk number 163: ch-gWidgetsIntro.Rnw:54-57
###################################################
options(prompt="> ")
options(continue="+ ")
options(width=80)


