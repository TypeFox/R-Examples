rzvvmd <-
  setRefClass("RzVVManageData",
              fields = c("main", "vvcore", "data"),
              methods = list(
                initialize            = function(...) {
                  initFields(...)
                  
                  # Data Management
                  rzVVSelectCases <- new("RzVVSelectCases", vvcore=vvcore, data=data)
                  rzVVDuplicateData <- new("RzVVDuplicateData", vvcore=vvcore, data=data)
                  notebook <- gtkNotebookNew()
                  notebook["tab-pos"] <- GtkPositionType["left"]
                  notebook$appendPage(rzVVSelectCases$getMain(), gtkLabelNew(gettext("Select Cases")))
                  notebook$appendPage(rzVVDuplicateData$getMain(), gtkLabelNew(gettext("Duplicate Dataset")))
                  
                  main <<- gtkWindowNew(show=FALSE)
                  main$setSizeRequest(400, 300)
                  main$setPosition(GtkWindowPosition["center-on-parent"])
                  main$setTitle(gettext("Manage Data"))
                  main$add(notebook)
                  
                  gSignalConnect(main, "delete-event", function(obj, ...) {
                    obj$hide()
                    return(TRUE)
                  })
                  
                },
                
                show = function(){
                  main$setTransientFor(vvcore$getWidget()$getToplevel())
                  main$showAll()
                },
                
                hide = function(){
                  main$hide()
                }
                
              )
  )
rzvvmd$accessors("main")
