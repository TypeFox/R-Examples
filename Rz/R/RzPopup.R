popup <- 
setRefClass("RzPopup",
  fields = c("uimanager", "context.menu", "a.recode"),
  methods = list(
    initialize            = function(...) {
      initFields(...)
      path <- file.path(path.package("Rz"), "ui", "popup.ui")
      
      a.recode    <<- gtkActionNew("Recode", gettext("Recode"))
      
      action.group  <- gtkActionGroupNew()
      action.group$setTranslationDomain("pkg-RGtk2")
      action.group$addAction(a.recode)

      uimanager     <<- gtkUIManagerNew()
      uimanager$insertActionGroup(action.group, 0)
      uimanager$addUiFromFile(path)
      context.menu <<- uimanager$getWidget("/PopupVV")
      
    },
    
    popupmenu = function(obj, event) {
        button <- event$GetButton()
        if (button==3) context.menu$popup(button=button,activate.time=event$getTime())
        return(TRUE)
    }
))
popup$accessors(c("uimanager", "context.menu", "a.recode"))
