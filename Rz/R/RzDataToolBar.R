datatoolbar <- 
  setRefClass("RzDataToolBar",
              fields = c("toolBar", "a.remove", "a.ch.name", "a.reload", "a.mngdata", "a.dataview"),
              methods = list(
                initialize            = function(...) {
                  initFields(...)
                  path1 <- file.path(path.package("Rz"), "ui", "datatoolbar.ui")
                  path2 <- file.path(path.package("Rz"), "images", "oxygen")
                  
                  a.remove    <<- gtkActionNew("RemoveDataSet", gettext("Remove Dataset"), gettext("Remove Dataset"))
                  a.ch.name   <<- gtkActionNew("ChageDataSetName", gettext("Change the Dataset Name"), gettext("Change the Dataset Name"))
                  a.reload    <<- gtkActionNew("ReloadFromGlobalEnv", gettext("Reload from Grobal Environment"), gettext("Reload from Grobal Environment"))
                  a.mngdata   <<- gtkActionNew("ManageData", gettext("Manage Data"), gettext("Manage Data"))
                  a.dataview  <<- gtkActionNew("DataView", gettext("View Data"), gettext("View Data"))
                  
                  a.remove$setIconFromFile(path2, "table_delete.png")
                  a.ch.name$setIconFromFile(path2, "edit-rename.png")
                  a.reload$setIconFromFile(path2, "table_refresh.png")
                  a.mngdata$setIconFromFile(path2, "manage_data.png")
                  a.dataview$setIconFromFile(path2, "table.png")
                  
                  action.group <- gtkActionGroupNew()
                  action.group$addAction(a.remove)
                  action.group$addAction(a.ch.name)
                  action.group$addAction(a.reload)
                  action.group$addAction(a.mngdata)
                  action.group$addAction(a.dataview)
                  
                  uimanager <- gtkUIManagerNew()
                  uimanager$insertActionGroup(action.group, 0)
                  uimanager$addUiFromFile(path1)
                  
                  toolBar <<- uimanager$getWidget("/DataToolBar")
                  toolBar["toolbar-style"] <<- GtkToolbarStyle["icons"]
                  toolBar$setIconSize(GtkIconSize["large-toolbar"])
                  nitems <- toolBar$getNItems()
                  for(i in seq_len(nitems)){
                    item <- toolBar$getNthItem(i-1)
                    child <- item$getChild()
                    child["can-focus"] <- FALSE
                  }
                  
                }
              ))
datatoolbar$accessors(c("toolBar", "a.remove", "a.ch.name", "a.reload", "a.mngdata", "a.dataview"))
