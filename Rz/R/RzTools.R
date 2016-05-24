tools <-
setRefClass("RzTools",
  fields = c("main", "variableView", "infoBar", "window", "datasetName", "dataHandler",
             "themeEditor", "tmpObj", "spinner", "RzPath", "initialized"),
  methods = list(
    initialize = function(...) {
      initFields(...)
      main         <<- NULL
      variableView <<- NULL
      infoBar      <<- NULL
      window       <<- NULL
      datasetName  <<- NULL
      themeEditor  <<- NULL
      tmpObj       <<- NULL
      dataHandler  <<- NULL
      spinner      <<- NULL
      initialized  <<- FALSE
    },
    
    clean = function(){
      main         <<- NULL
      variableView <<- NULL
      infoBar      <<- NULL
      window       <<- NULL
      datasetName  <<- NULL
      dataHandler  <<- NULL
      spinner      <<- NULL
    },
    
    runDialog = function(message, type = c("info", "warning", "question", "error")){
      type <- match.arg(type)
      buttons.type <- c(info="ok", warning="close", question="ok-cancel", error="close")
      dialog <- gtkMessageDialogNew(window, "destroy-with-parent",
                                    GtkMessageType[type], GtkButtonsType[buttons.type[type]],
                                    message)
      response <- dialog$run()
      dialog$hide()        
      return(response)
    },
    
    addData = function(data.set, name){
      if (is.data.frame(data.set)) {
        data.set2 <- data.set(data.set)
        names(data.set2) <- colnames(data.set)
        data.set <- data.set2
      }
      
      if (all(class(data.set) == "data.set")) {
        data <- rzdata$new(file.path=NULL, data.set=data.set, original.name=name)
        if (is.null(main)) {
          data.collection.obj$addData(data)
        } else {
          main$addData(data)
        }
      } else {
        stop("invalid 'type' (", class(data.set)[1], ") of argument")
      }
    },
    
    addItem = function(item, name = as.character(substitute(item)), data.set.name = NULL, description = name,
                       measurement = c("auto", "nominal", "ordinal", "interval", "ratio"),
                       overwrite = FALSE, ask = FALSE) {
      if (is.null(main)) {
        stop("Please start Rz.")        
      } else {
        main$addItem(item = item, name = name, data.set.name, description = description,
                     measurement = measurement, overwrite = overwrite, ask = ask)
      }
            
    },
    
    reloadData = function(data.set.name=NULL, ask = TRUE) {
      if (is.null(main)) {
        stop("Please start Rz.")
      } else {
        main$reloadData(data.set.name, ask = ask)
      }
    },
    
    sync = function(name, obj){
      tmpObj <<- obj
      eval(parse(text=paste(name, " <- Rz:::rzTools$getTmpObj()")), envir=.GlobalEnv)
    }
  )
)
tools$accessors("main", "variableView", "infoBar", "window", "datasetName",
                "themeEditor", "tmpObj", "dataHandler", "spinner", "RzPath", "initialized")

