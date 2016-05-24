vartoolbar <- 
  setRefClass("RzVarToolBar",
              fields = c("toolBar", "a.delete", "a.duplicate", "a.select", "a.unselect", "a.onlyselec", "signal.id"),
              methods = list(
                initialize            = function(...) {
                  initFields(...)
                  signal.id <<- list()
                  path1 <- file.path(path.package("Rz"), "ui", "vartoolbar.ui")
                  path2 <- file.path(path.package("Rz"), "images", "oxygen")
                  
                  a.delete    <<- gtkActionNew("Delete", gettext("Delete Variables"), gettext("Delete Variables"))
                  a.duplicate <<- gtkActionNew("Duplicate", gettext("Duplicate Variables"), gettext("Duplicate Variables"))
                  a.select    <<- gtkActionNew("SelectAll", gettext("Select All Variables"), gettext("Select All Variables"))
                  a.unselect  <<- gtkActionNew("UnselectAll", gettext("Unselect All Variables"), gettext("Unselect All Variables"))
                  a.onlyselec <<- gtkToggleActionNew("ToggleSelectedView", gettext("Show Selected Variables Only"), gettext("Show Selected Variables Only"))
                  
                  a.delete$setIconFromFile(path2, "table_row_delete.png")
                  a.duplicate$setIconFromFile(path2, "table_row_insert.png")
                  a.select$setIconFromFile(path2, "checkbox.png")
                  a.unselect$setIconFromFile(path2, "nocheck.png")
                  a.onlyselec$setIconFromFile(path2, "view-filter.png")
                  
                  action.group  <- gtkActionGroupNew()
                  action.group$addAction(a.delete)
                  action.group$addAction(a.duplicate)
                  action.group$addAction(a.select)
                  action.group$addAction(a.unselect)
                  action.group$addAction(a.onlyselec)
                  
                  uimanager     <- gtkUIManagerNew()
                  uimanager$insertActionGroup(action.group, 0)
                  uimanager$addUiFromFile(path1)
                  
                  toolBar     <<- uimanager$getWidget("/VarToolBar")
                  toolBar$setShowArrow(TRUE)
                  toolBar["toolbar-style"] <<- GtkToolbarStyle["icons"]
                  toolBar$setIconSize(GtkIconSize["large-toolbar"])
                  nitems <- toolBar$getNItems()
                  for(i in seq_len(nitems)){
                    item <- toolBar$getNthItem(i-1)
                    child <- item$getChild()
                    child["can-focus"] <- FALSE
                  }
                  toolBar$showAll()
                  
                },
                
                setSignalId = function(name, id) {
                  signal.id[[name]] <<- id
                },
                
                setFiltered = function(active) {
                  gSignalHandlerBlock(a.onlyselec, signal.id[["onlyselec"]])
                  a.onlyselec$setActive(active)
                  gSignalHandlerUnblock(a.onlyselec, signal.id[["onlyselec"]])
                }
              ))
vartoolbar$accessors(c("toolBar", "a.delete", "a.duplicate", "a.select", "a.unselect", "a.onlyselec"))
