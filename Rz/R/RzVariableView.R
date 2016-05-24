variable.view <- 
setRefClass("RzVariableView",
  fields = c("data", "win", "main", "liststore", "summaries", "summaries.subset",
             "nominalpix", "ordinalpix", "intervalpix", "ratiopix", "vvs"),
  methods = list(
    initialize  = function(...) {
      initFields(...)
      nominalpix    <<- gdkPixbufNewFromFile(file.path(path.package("Rz"), "images/oxygen/cat.png"     ))$retval
      ordinalpix    <<- gdkPixbufNewFromFile(file.path(path.package("Rz"), "images/oxygen/order.png"   ))$retval
      intervalpix   <<- gdkPixbufNewFromFile(file.path(path.package("Rz"), "images/oxygen/interval.png"))$retval
      ratiopix      <<- gdkPixbufNewFromFile(file.path(path.package("Rz"), "images/oxygen/ratio.png"   ))$retval
      
      liststore <<- gtkListStoreNew("character", "logical", "character", "character", "character", "GdkPixbuf", "character", "character")
      vars      <-  data$getVariableNames()
      var.labs  <-  data$getVariableLabels()
      msr       <-  data$getMeasurement()
      val.labs  <-  data$getValueLabels()
      miss.val  <-  data$getMissingValues()
      summaries <<- data$getSummaries()
      summaries.subset <<- NULL
      for ( i in seq_len(data$ncol()) ) {
        iter <- liststore$append()$iter
        liststore$set(iter,
                      column.definition["index"], i,
                      column.definition["select"], FALSE,
                      column.definition["vars"], vars[i],
                      column.definition["var.labs"], var.labs[i],
                      column.definition["msr"], msr[i],
                      column.definition["msr.image"], .self$msrPix(msr[i]),
                      column.definition["val.labs"], val.labs[i],
                      column.definition["missing"], miss.val[i])
      }
      
      rzVVMain         <- new("RzVVMain", data=data, liststore=liststore, vv=.self)
      adj <- rzVVMain$getAdj()
      rzVVDescriptives <- new("RzVVDescriptives", data=data, liststore=liststore, vv=.self, adj=adj)
      rzVVPlot         <- new("RzVVPlot", data=data, liststore=liststore, vv=.self, adj=adj)
      rzVVManageData   <- new("RzVVManageData", vvcore=rzVVMain, data=data)
      vvs <<- list(rzVVMain, rzVVDescriptives, rzVVPlot)
      
      # notebook
      hbox.datatoolbar <- gtkHBoxNew(spacing=2)
      rzDataToolBar <- new("RzDataToolBar")
      toolbar <- rzDataToolBar$getToolBar()
      toolbar$showAll()
      toolbar$setShowArrow(FALSE)
      hbox.datatoolbar$packEnd(toolbar, expand=TRUE)
      
      main <<- gtkNotebookNew()
      if(is.null(gtkCheckVersion(2, 20, 0))) {
        main$setActionWidget(hbox.datatoolbar, GtkPackType["end"])
      }
      
      main$appendPage(rzVVMain$getMain(), gtkLabelNew(gettext("Variables")))
      main$appendPage(rzVVDescriptives$getMain(), gtkLabelNew(gettext("Descriptive statistics")))
      main$appendPage(rzVVPlot$getMain(), gtkLabelNew(gettext("Plot")))
      
      gSignalConnect(main, "switch-page", function(obj, page, page_num){
        
      })
      gSignalConnect(rzDataToolBar$getA.ch.name(), "activate", .self$onChangeDataSetName)
      gSignalConnect(rzDataToolBar$getA.remove(), "activate", .self$onRemove)
      gSignalConnect(rzDataToolBar$getA.reload(), "activate", .self$onReload)
      gSignalConnect(rzDataToolBar$getA.mngdata(), "activate", function(...){
        rzVVManageData$show()
      })
      gSignalConnect(rzDataToolBar$getA.dataview(), "activate", function(...){
        new("RzDataView", RzData=data)
      })
      gSignalConnect(main, "hide", function(...){
        rzVVManageData$hide()
      })
      
    },
    
    setFiltered = function(action){
      if(action$getActive()) {
        lapply(vvs, function(x) x$setFiltered(TRUE))
      } else {
        lapply(vvs, function(x) x$setFiltered(FALSE))
      }
    },
    
    msrPix = function(msr){
      pix <- NULL
      if(msr=="nominal"){
        pix <- nominalpix
      } else if(msr=="ordinal"){
        pix <- ordinalpix
      } else if(msr=="interval"){
        pix <- intervalpix
      } else if(msr=="ratio"){
        pix <- ratiopix
      }
      return(pix)
    },
    
    showAll = function(){
      .self$changeFont()
      main$showAll()
    },

    hide = function(){
      main$hide()
    },
    
    reload = function(){
      on.exit(spinStop())
      if (data$getSubset.on()) {
        if (main$getRealized()) {
          dialog <- gtkMessageDialogNew(rzTools$getWindow(), "destroy-with-parent",
                                        "error", "close", gettext("Cannot reload while Select Cases enabled."))
          dialog$run()
          dialog$hide()
          return()
          
        } else {
          stop("Cannot reload while \"Select Cases\" enabled.")          
        }
      }
      spinStart()
      result <- data$reloadFromGlobalEnv()
      
      if (result != TRUE) {
        if (main$getRealized()) {
          dialog2 <- gtkMessageDialogNew(win, "destroy-with-parent",
                                         GtkMessageType["error"], GtkButtonsType["close"],
                                         gettextf("\"%s\" isn't a data.set or doesn't exist.", result))
          dialog2$run()
          dialog2$hide()
        } else {
          stop("\"", result, "\" isn't a data.set or doesn't exist.")          
        }
      }
      
      iter <- liststore$getIterFirst()
      selects <- logical(0)
      while(iter$retval){
        select  <- liststore$getValue(iter$iter, column.definition["select"])$value
        selects <- c(selects, select)
        iter$retval <- liststore$iterNext(iter$iter)
      }
      diff <- data$ncol() - length(selects)
      if(diff > 0){
        selects <- c(selects, rep(FALSE, diff))        
      } else if(diff < 0){
        selects <- rep(FALSE, data$ncol())
      }
      vars      <-  data$getVariableNames()
      var.labs  <-  data$getVariableLabels()
      msr       <-  data$getMeasurement()
      val.labs  <-  data$getValueLabels()
      miss.val  <-  data$getMissingValues()
      summaries <<- data$getSummaries()
      liststore$clear()
      for ( i in seq_len(data$ncol()) ) {
        iter <- liststore$append()$iter
        liststore$set(iter,
                      column.definition["index"], i,
                      column.definition["select"], selects[i],
                      column.definition["vars"], vars[i],
                      column.definition["var.labs"], var.labs[i],
                      column.definition["msr"], msr[i],
                      column.definition["msr.image"], .self$msrPix(msr[i]),
                      column.definition["val.labs"], val.labs[i],
                      column.definition["missing"], miss.val[i])
      }
    },
    
    # actions
    onChangeDataSetName = function(...){
      rzTools$getMain()$changeDataSetName()
    },
    
    onRemove = function(...){
      rzTools$getMain()$remove()
    },
    
    onReload = function(action){
      dialog <- gtkMessageDialogNew(main$getToplevel(), "destroy-with-parent",
                                    GtkMessageType["question"], GtkButtonsType["ok-cancel"],
                                    gettext("Are you sure you want to do that?"))
      response <- dialog$run()
      dialog$hide()
      
      if(response==GtkResponseType["ok"]){
        .self$reload()
      }
    },
    
    getView = function() return(main),
    
    setSubsetSummaries = function(){
      summaries.subset <<- data$getSummaries(subset=TRUE)
    },
    
    setSummary = function(row){
      summaries[row] <<- data$getSummary(row)
      if (data$getSubset.on() & nzchar(data$getSubset.condition()))
        summaries.subset[row] <<- data$getSummary(row, subset=TRUE)
    },
    
    changeFont = function(){
      lapply(vvs, function(x) x$changeFont())
    }
    
  )
)
variable.view$accessors(c("liststore", "data", "summaries", "summaries.subset"))
