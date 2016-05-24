main <-
setRefClass("RzMain",
  fields = c("win",
             "main.hpaned", "main.vpaned", "main.view", "variable.view",
             "info.bar", "message.label",
             "variable.view.list",
             "rzMenu", "rzDataHandler", "rzDataSetIO",
             "view.box"),
  methods = list(
    initialize            = function(...) {
      initFields(...)
      if (! exists("theme_rz", envir=.GlobalEnv)) {
        rzTools$sync("theme_rz", theme_grey)        
      }
      
      if (rzSettings$getAutosave() && ! rzTools$getInitialized()) loadSession()
      rzTools$setInitialized(TRUE)
      

      settings <- gtkSettingsGetDefault()
      settings$setStringProperty("gtk-font-name", rzSettings$getGlobalFont(), NULL)
      gtkRcReparseAll()
      rzDataSetIO    <<- new("RzDataSetIO")
      win            <<- gtkWindowNew(show=FALSE)
      rzTools$setWindow(win)
      info.bar      <<- gtkInfoBarRzNew()
      message.label <<- gtkLabelNew()
      info.bar$addButton("gtk-ok", GtkResponseType["ok"])
      info.bar$setMessageType(GtkMessageType["info"])
      info.bar$getContentArea()$packStart(message.label, expand=FALSE, fill=FALSE)
      info.bar$hide()
      info.bar$setNoShowAll(TRUE)
      rzTools$setInfoBar(info.bar)
      spinner       <- gtkSpinnerNew(show=FALSE)
      spinner$setSizeRequest(22, 22)
      rzTools$setSpinner(spinner)
      
      main.hpaned    <<- gtkHPanedNew()
      main.vpaned    <<- gtkVPanedNew()
      main.view      <<- gtkVBoxNew()
      
      main.vpaned$pack1(main.view, resize=TRUE)
      main.vpaned$setPosition(280)
      
      main.hpaned$pack1(main.vpaned, resize=TRUE)
      main.hpaned$setPosition(450)
      
      variable.view <<- NULL
      variable.view.list <<- list()
      
      win["title"] <<- "Rz"
      win$setDefaultSize(800, 700)
      win["window-position"] <<- GtkWindowPosition["center"]
      
      gSignalConnect(info.bar, "response", .self$onInfoBarResponsed)
      
      rzMenu <<- new("RzMenu")
      accel.group <- rzMenu$getUimanager()$getAccelGroup()
      
      win$addAccelGroup(accel.group)
      
      gSignalConnect(rzMenu$getA.open(),      "activate", .self$onOpen, "rzdata")
      gSignalConnect(rzMenu$getA.save(),      "activate", .self$onSave)
      gSignalConnect(rzMenu$getA.importsav(), "activate", .self$onOpen, "spss")
      gSignalConnect(rzMenu$getA.importpor(), "activate", .self$onOpen, "spss.por")
      gSignalConnect(rzMenu$getA.importdta(), "activate", .self$onOpen, "stata")
      gSignalConnect(rzMenu$getA.importt(),   "activate", .self$onOpen, "text")
      gSignalConnect(rzMenu$getA.importg(),   "activate", .self$onImportFromGlobalEnv)
      gSignalConnect(rzMenu$getA.quit(),      "activate", win$destroy)
      gSignalConnect(rzMenu$getA.settings(),  "activate", .self$onSetting)
      gSignalConnect(rzMenu$getA.tutorial(),  "activate", function(...) browseURL(gettext("http://m884.jp/RzTutorial.html")))
      gSignalConnect(rzMenu$getA.load.sample(), "activate", .self$onLoadSample)
      
      rzDataHandler <<- new("RzDataHandler")
      rzTools$setDataHandler(rzDataHandler)
      gSignalConnect(rzDataHandler$getData.set.list.combo(), "changed", .self$onDatasetChanged)
      
      vbox      <- gtkVBoxNew()
      hbox.menu <- gtkHBoxNew()
      hbox.menu$packStart(rzMenu$getMenu.bar(), expand=TRUE, fill=TRUE)
      hbox.menu$packEnd(spinner, expand=FALSE, padding=5)
      hbox.menu$packEnd(rzDataHandler$getData.set.list.combo(), expand=FALSE)
      
      vbox$packStart(hbox.menu,    expand=FALSE, fill=FALSE)
      vbox$packStart(info.bar,     expand=FALSE, fill=FALSE)
      vbox$packStart(main.hpaned,  expand=TRUE , fill=TRUE)
      win$add(vbox)
      win$show()
      win$present()

      gSignalConnect(win, "destroy", function(...){
        rzTools$clean()
        if(!is.null(rzSettings$getEmbededDeviceOn()) && rzSettings$getEmbededDeviceOn()){
          try(dev.off(), silent=TRUE)
        }
        if (rzSettings$getAutosave()) saveSession()
      })
      
      rzDataHandler$changeData(0)
      
    },
    
    show = function(){
      win$show()
      win$present()
    },
    
    # actions
    onSetting = function(action){
      rzSettings$runDialog(win)
      if(!is.null(variable.view))variable.view$changeFont()
      gtkRcReparseAll()
    },
    
    onOpen                = function(action, type){
      on.exit(spinStop())
      if (type=="text") {
        data <- rzDataSetIO$importFromText(win)
      } else {
        data <- rzDataSetIO$open(win, type, info.bar)
      }
      if(!is.null(data)) {
        rzDataHandler$addData(data)
      }
    },
    
    onSave                = function(action){
      on.exit(spinStop())
      data <- rzDataHandler$getCurrentData()
      rzDataSetIO$save(win, data)
    },
    
    onImportFromGlobalEnv = function(action){
      on.exit(spinStop())
      data <- rzDataSetIO$importFromGlobalEnv(win)
      if(!is.null(data)) {
        rzDataHandler$addData(data)
      }
    },
    
    onLoadSample = function(action){
      on.exit(spinStop())
      spinStart()
      nes1948.por <- UnZip("anes/NES1948.ZIP","NES1948.POR",package="memisc")
      nes1948 <- spss.portable.file(nes1948.por)
      sample.data.set <- as.data.set(nes1948)
      data <- new("RzData", file.path=NULL, data.set=sample.data.set,
                  original.name=gettext("NES1948 [Sample Dataset in memisc]"))
      rzDataHandler$addData(data)
    },
    
    changeDataSetName   = function(){
      onActivate <- function(entry, dialog){
        dialog$response(GtkResponseType["ok"])
      }
      onResponse <- function (dialog, response.id, entry) {
        if (response.id==GtkResponseType["ok"]) {
          new.name <- localize(entry$getText())
          result <- rzDataHandler$changeDataSetName(data.set.name, new.name)
          if(result$result){
            rm(list=data.set.name, envir=.GlobalEnv)
            rm(list=sprintf("%s.ds", data.set.name), envir=.GlobalEnv)
            match <- which(names(variable.view.list)==data.set.name)
            names(variable.view.list)[match] <<- new.name
            dialog$hide()            
          } else {
            dialog2 <- gtkMessageDialogNew(dialog, "destroy-with-parent",
                                           GtkMessageType["error"],
                                           GtkButtonsType["close"],
                                           result$message)
            dialog2$run()
            dialog2$hide()
            dialog$run()
            return()
          }
        } else {
          dialog$hide()
        }        
      }
      
      data.set.name <- rzDataHandler$getCurrentDataSetName()
      if ( length(data.set.name) == 0 ) return()
      txt <- paste(gettextf("Old Name"),":", sprintf(" <b>%s</b>", data.set.name), sep="")
      label1 <- gtkLabelNew()
      label1$setMarkup(txt)
      label1$show()
      label2 <- gtkLabelNew(gettext("New Name"))
      entry1 <- gtkEntryNew()
      hbox1 <- gtkHBoxNew()
      hbox2 <- gtkHBoxNew()
      hbox1$packStart(label1, expand = FALSE, padding = 5)
      hbox2$packStart(label2, expand = FALSE, padding = 5)
      hbox2$packStart(entry1, expand = TRUE, fill = TRUE, padding = 5)
      
      vbox1 <- gtkVBoxNew()
      vbox1$packStart(hbox1, expand = FALSE, padding = 2)
      vbox1$packStart(hbox2, expand = FALSE, padding = 2)
      
      dialog <- gtkDialogNewWithButtons(title=gettext("Change the Name of Current Dataset"), parent=win, flags=c("modal", "destroy-with-parent"),
                                        "gtk-ok", GtkResponseType["ok"], 
                                        "gtk-cancel", GtkResponseType["cancel"],
                                        show=FALSE)
      dialog$setResizable(FALSE)
      dialog$setDefaultSize(200, 130)
      dialog[["vbox"]]$packStart(vbox1, expand = FALSE, padding = 2)
      gSignalConnect(entry1, "activate", onActivate, dialog)
      gSignalConnect(dialog, "response", onResponse, entry1)
      response <- dialog$run()
    },
    
    remove              = function(){
      data.set.name <- rzDataHandler$getCurrentDataSetName()
      if ( length(data.set.name) == 0 ) return()
      dialog <- gtkMessageDialogNew(win, "destroy-with-parent",
                                    GtkMessageType["question"], GtkButtonsType["ok-cancel"],
                                    gettext("Are you sure you want to do that?"))
      response <- dialog$run()
      dialog$hide()
      
      if(response==GtkResponseType["ok"]){
        ind <- rzDataHandler$removeCurrentData()
        variable.view$getView()$destroy()
        variable.view <<- NULL
        variable.view.list[data.set.name] <<- NULL
        rzTools$setVariableView(variable.view)
        rzDataHandler$changeData(ifelse(ind==0, 0, ind-1))
      }
    },
    
    onRevert = function(action){
      on.exit(spinStop())
      spinStart()
      if(!is.null(variable.view)){
        dialog <- gtkMessageDialogNew(win, "destroy-with-parent",
                                      GtkMessageType["question"], GtkButtonsType["ok-cancel"],
                                      gettext("Are you sure you want to do that?"))
        response <- dialog$run()
        dialog$hide()
        
        if(response==GtkResponseType["ok"]){
          data <- rzDataHandler$getCurrentData()
          data$revert()
          variable.view$reload()
        }
      }
    },
    
    onDatasetChanged      = function(combo){
      on.exit(spinStop())
      spinStart()
      data.set.name <- rzDataHandler$getCurrentDataSetName()
      if ( length(data.set.name) == 0 ) return()
      
      if ( !is.null(variable.view) ) {
        variable.view$hide()
      }
      
      if ( ! data.set.name %in% names(variable.view.list) ) {
        variable.view <<- new("RzVariableView", data=rzDataHandler$getData(data.set.name),
                              win=win)
        variable.view.list[[data.set.name]] <<- variable.view
        main.view$packStart(variable.view$getView())
      } else {
        variable.view <<- variable.view.list[[data.set.name]]
        variable.view$showAll()
      }
      rzTools$setVariableView(variable.view)
      rzDataHandler$sync(data.set.name)
      
      win["title"] <<- paste("Rz -", data.set.name)
    },
    
    onInfoBarResponsed    = function(widget, response.id){
      info.bar$hide()
      
    },
    
    # scripting interface
    addData = function(data){
      rzDataHandler$addData(data)
    },
    
    addItem = function(item, name = as.character(substitute(item)), data.set.name = NULL, description = name,
                       measurement = c("auto", "nominal", "ordinal", "interval", "ratio"),
                       overwrite = FALSE, ask = FALSE) {
      measurement <- match.arg(measurement)
      item <- as.item(item)
      if (!is.null(description))
        description(item) <- description
      if (measurement != "auto")
        measurement(item) <- measurement
      
      vv.tmp <- NULL
      if (is.null(data.set.name)) {
        vv.tmp <- variable.view
      } else {
        vv.tmp <- variable.view.list[[data.set.name]]
      }
      if (! is.null(vv.tmp)) {
        response <- 1
        if (ask) {
          response <- menu(c(gettext("yes"), gettext("no")),
                           title = gettext("Are you sure you want to do that?"))
        }
        
        if (response == 1) {
          data.tmp <- vv.tmp$getData()
          exist <- data.tmp$existVar(name)
          if (exist && !overwrite) {
            response <- menu(c(gettext("yes"), gettext("no")),
                             title = gettextf("\"%s\" already exists. Are you sure you want to overwrite?",
                                              name))
            if (response != 1)
              return()
          }
          data.tmp$addItem(item, name)
          vv.tmp$reload()
        }
        
      } else {
        if (is.null(data.set.name)) {
          stop("Please select a dataset on Rz or specify \"data.set.name\".")
        } else {
          stop("The dataset named \"", data.set.name, "\" doesn't exist.")
        }
      }
    },
    
    reloadData = function(data.set.name=NULL, ask = TRUE){
      vv.tmp <- NULL
      if (is.null(data.set.name)) {
        vv.tmp <- variable.view
      } else {
        vv.tmp <- variable.view.list[[data.set.name]]
      }
      if (! is.null(vv.tmp)) {
        response <- 1
        if (ask) {
          response <- menu(c(gettext("yes"), gettext("no")),
                           title = gettext("Are you sure you want to do that?"))
        }
        
        if (response == 1) {
          vv.tmp$reload()
        }
        
      } else {
        if (is.null(data.set.name)) {
          stop("Please select a dataset on Rz or specify \"data.set.name\".")
        } else {
          stop("The dataset named \"", data.set.name, "\" doesn't exist.")
        }
      }
    }
  )
)
main$accessors("variable.view")
