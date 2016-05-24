
###################################################
### code chunk number 122: Dialogs.Rnw:121-147
###################################################
window <- tktoplevel()
tkwm.title(window, "File menu example")
menu_bar <- tkmenu(window)
tkconfigure(window, menu = menu_bar)
file_menu <- tkmenu(menu_bar)
tkadd(menu_bar, "cascade", label="File", menu = file_menu)

tkadd(file_menu,"command", label = "Source file...",
      command =  function() {
        file_name <- tkgetOpenFile(filetypes=
                        "{{R files} {.R}} {{All files} *}")
        if(file.exists(file_name <- as.character(file_name)))
           source(tclvalue(file_name))
      })
tkadd(file_menu, "command", label = "Save workspace as...",
      command = function() {
        file_name <- tkgetSaveFile(defaultextension = "Rsave")
        if(nchar(fname <- as.character(file_name)))
          save.image(file = file_name)
      })
tkadd(file_menu, "command", label="Set working directory...",
      command = function() {
        dir_name <- tkchooseDirectory()
        if(nchar(dir_name <- as.character(dir_name)))
          setwd(dir_name)
      })
