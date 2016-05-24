
###################################################
### code chunk number 181: Dialogs.Rnw:487-494
###################################################
name_filter <- paste("R files (*.R .RData)",
                    "Sweave files (*.Rnw)",
                    "All files (*.*)", 
                    sep=";;")
##
filenames <- Qt$QFileDialog$getOpenFileNames(NULL, 
             "Open file(s)...", getwd(), name_filter)


print(filenames)
