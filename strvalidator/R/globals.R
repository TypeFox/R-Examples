# Hack to avoid NOTES in R CMD check
# Hadley does not seem to like it: 
# http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
# Alternatively '::' can be used: https://stat.ethz.ch/pipermail/r-devel/2007-June/046048.html

## Needed to avoid notes when using gWidgets.
if (base::getRversion() >= "2.15.1") {
  utils::globalVariables(c("gwindow","ggroup","glayout","glabel","gdroplist",
                           "addHandlerChanged","svalue","svalue<-","gmessage",
                           "gframe","gedit","gbutton","enabled<-","dispose",
                           "gradio","gcheckbox","visible","visible<-",
                           "gtable","gspinbutton","gtext","delete",
                           "gfilebrowse","gbasicdialog","gexpandgroup",
                           "addHandlerDoubleclick","addDropSource","addDropTarget",
                           "blockHandler","addHandlerFocus","ginput","gfile",
                           "gnotebook"))
}

## Needed to avoid notes when using data.table in functions:
## calculateAllele, calculateAT, calculateAT_gui, calculateAT6,
## calculateAT6_gui, plotBalance_gui, calculateLb, calculateRatio,
## generateEPG
if (base::getRversion() >= "2.15.1") {
  utils::globalVariables(c("Height","Blocked","Sample.File.Name","Dye",
                           "Hb","Marker","Lb", "H", "Amount", "Size",
                           "Sample.Name", "TPH", "MPH", "TPPH", "MTPH",
                           "Group", "Allele", "Id"))
}