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


