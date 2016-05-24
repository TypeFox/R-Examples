rzplot <- 
setRefClass("RzPlot",
  fields = c("info.bar", "main", "notebook", "data", "model", "constructed", "vvcore",
             "geom", "stratum", "misc", "save",
             "button.prev", "button.next",
             "p.current", "p.list", "p.current.num", "parent", "parent.win", "rzPlotScript"),
  methods = list(
    initialize  = function(...) {
      initFields(...)
      constructed   <<- FALSE
      info.bar      <<- rzTools$getInfoBar()
      parent.win    <<- rzTools$getWindow()
      parent        <<- NULL
      p.current     <<- NULL
      model         <<- NULL
      data          <<- NULL
      rzPlotScript  <<- new("RzPlotScript")
      p.list        <<- list()
      p.current.num <<- 1
      save          <<- new("RzPlotSave")
      
      
      notebook <<- gtkNotebookNew()
      button.detach <- gtkButtonNewWithLabel(gettext("Detach"))
      notebook$setTabPos(GtkPositionType["top"])
      notebook$setScrollable(TRUE)

      button.save          <- save$getButton.save()
      button.execute       <- gtkButtonNewFromStock(GTK_STOCK_EXECUTE)
      
      image  <- gtkImageNewFromStock(GTK_STOCK_CLEAR, GtkIconSize["button"])
      button.clear          <- gtkButtonNew()
      button.clear["sensitive"] <- TRUE
      button.clear$setImage(image)
      
      image  <- gtkImageNewFromStock(GTK_STOCK_GO_BACK, GtkIconSize["button"])
      button.prev          <<- gtkButtonNew()
      button.prev["sensitive"] <<- FALSE
      button.prev$setImage(image)
      
      image  <- gtkImageNewFromStock(GTK_STOCK_GO_FORWARD, GtkIconSize["button"])
      button.next          <<- gtkButtonNew()
      button.next["sensitive"] <<- FALSE
      button.next$setImage(image)
      
      button.box.history <- gtkHBoxNew()
      button.box.history$packStart(button.prev, expand=FALSE)
      button.box.history$packStart(button.next, expand=FALSE)
      
      entry.script <- gtkEntryNew()
      entry.script$setIconFromStock(GtkPositionType["left"], GTK_STOCK_JUSTIFY_LEFT)
      entry.script$setText("View Script")
      entry.script$setWidthChars(15)
      entry.script$setEditable(FALSE)
      
      gSignalConnect(entry.script, "motion-notify-event", function(...){
        popup <- gtkWindowNew(type=GtkWindowType["toplevel"], show=FALSE)
        popup$setDecorated(FALSE)
        popup$setTransientFor(parent.win)
        alloc <- entry.script$getAllocation()
        popup$setDefaultSize(alloc$allocation$width + 1, alloc$allocation$height + 1)
        
        table <- gtkTextTagTableNew()
        tag <- gtkTextTagNew("first-line")
        tag["pixels-above-lines"] <- 10
        table$add(tag)
        tag <- gtkTextTagNew("last-line")
        tag["pixels-below-lines"] <- 10
        table$add(tag)
        
        tw     <- gtkTextViewNew()
        tw$modifyFont(pangoFontDescriptionFromString(rzSettings$getMonospaceFont()))
        
        buffer <- gtkTextBufferNew(table)
        if(!is.null(data)) buffer$setText(rzPlotScript$getScript())
        first1 <- buffer$getStartIter()$iter
        first2 <- buffer$getStartIter()$iter
        first2$forwardToLineEnd()
        last1  <- buffer$getEndIter()$iter
        last2  <- buffer$getIterAtLine(-1)$iter
        buffer$applyTagByName("first-line", first1, first2)
        buffer$applyTagByName("last-line" , last1 , last2)
        tw$setBuffer(buffer)
        tw$setEditable(FALSE)
        tw$setLeftMargin(10)
        tw$setRightMargin(10)
        frame <- gtkFrameNew()
        frame$add(tw)
        frame$setShadowType(GtkShadowType["in"])
        popup$add(frame)
        
        
        geo <- parent.win$getWindow()$getGeometry()
        relative.coord <- entry.script$translateCoordinates(parent.win, 0, 0)
        geo$x <- geo$x + relative.coord$dest.x
        geo$y <- geo$y + relative.coord$dest.y
        popup$move(geo$x, geo$y)
        popup$showAll()
        tw$grabFocus()
        
        gSignalConnect(tw, "focus-out-event", function(widget, event){popup$destroy(); return(FALSE)})
        gSignalConnect(popup, "leave-notify-event", function(widget, event){
          size <- popup$getSize()
          cond <- c(event$x, event$y ,size$width - event$x, size$height - event$y)
          if (any(cond < 0))  popup$destroy()
          return(FALSE)
        })
        return(FALSE)
      })
      
      button.box1 <- gtkHBoxNew(spacing=5)
      button.box1$packStart(entry.script    , expand=FALSE)
      button.box1$packEnd(button.save       , expand=FALSE)
      button.box1$packEnd(button.execute    , expand=FALSE)
      button.box1$packEnd(button.box.history, expand=FALSE)
      button.box1$packEnd(button.clear      , expand=FALSE)
      
      vbox <- gtkVBoxNew()
      vbox$packStart(button.box1, expand=FALSE, padding=5)
      if(is.null(gtkCheckVersion(2, 20, 0))) {
#        notebook$setActionWidget(button.detach, GtkPackType["end"])
      }      
      vbox$packStart(notebook, expand=TRUE, fill=TRUE)
      
      main <<- gtkVPanedNew()
      main$pack2(vbox, resize=FALSE)
      
      gSignalConnect(button.detach, "clicked", .self$onDetach)
      
      gSignalConnect(button.execute, "clicked", function(button){
        if(!is.null(data)) .self$onPlot()
      })
      
      gSignalConnect(button.clear, "clicked", function(button){
        geom$clear()
        stratum$clear()
        misc$clear()
      })
      
      gSignalConnect(button.prev, "clicked", function(button){
        p.current.num <<- p.current.num + 1
        p.current <<- p.list[[p.current.num]]
        suppressWarnings(print(p.current))
        button.next["sensitive"] <<- TRUE
        if(p.current.num==length(p.list)) button["sensitive"] <- FALSE        
      })
      
      gSignalConnect(button.next, "clicked", function(button){
        p.current.num <<- p.current.num - 1
        p.current <<- p.list[[p.current.num]]
        suppressWarnings(print(p.current))
        button.prev["sensitive"] <<- TRUE
        if(p.current.num==1) button["sensitive"] <- FALSE        
      })
      
    },
    
    construct  = function(...) {
      
      main$hide()
      
      geom     <<- new("RzPlotGeom", rzPlotScript=rzPlotScript)
      stratum  <<- new("RzPlotStratum", rzPlotScript=rzPlotScript)
      misc     <<- new("RzPlotMisc", rzPlotScript=rzPlotScript)
      
      geom$completionSetModel(model)
      stratum$completionSetModel(model)
      
      # container
      notebook$foreach(gtkWidgetDestroy)
      notebook$appendPage(geom$getMain(), tab.label=gtkLabelNew(gettext("Main")))
      notebook$appendPage(stratum$getStratumPage(), tab.label=gtkLabelNew(gettext("Sub-group")))
      notebook$appendPage(stratum$getFacetPage(), tab.label=gtkLabelNew(gettext("Facet")))
      notebook$appendPage(misc$getAxisPage(), tab.label=gtkLabelNew(gettext("Axis")))
      notebook$appendPage(misc$getThemePage(), tab.label=gtkLabelNew(gettext("Theme")))
            
      if(rzSettings$getUseEmbededDevice()){
        if(require(cairoDevice)) {
          try(dev.off(), silent=TRUE)
          plot.area <- gtkDrawingArea()
          main$pack1(plot.area, resize=TRUE)          
          asCairoDevice(plot.area)
          rzSettings$setEmbededDeviceOn(TRUE)
          main$setPosition(300)
        } else {
          cat("cairoDevice package isn't installed", fill=TRUE)
          rzSettings$setEmbededDeviceOn(FALSE)          
        }
      } else {
        rzSettings$setEmbededDeviceOn(FALSE)
      }
      
      main$showAll()
      if (!constructed & rzSettings$getEmbededDeviceOn()) {
        p <- ggplot()
        p <- p + annotate("text", label = "Plot Area", x=10, y=10, size=18)
        p <- p + theme(axis.ticks=element_blank(), axis.text=element_blank(),
                       axis.title=element_blank())
        print(p)
      }
      
      gSignalConnect(save$getDialog(), "response", .self$onSave)
      
      constructed <<- TRUE
    },
    
    onSave = function(dialog, response.id){
      save$onSave(response.id, p.current, data)
    },
    
    setX = function(col){
      variableNames <- data$getVariableNames()
      geom$setX(variableNames[col])
    },
    
    onPlot = function(){      
      script <- rzPlotScript$getScript()
      env <- new.env()
      env[[data$getDataSetName()]] <- data$getData.frame()

      con <- textConnection("str", open="w", local=TRUE)
      str <- ""
#      cat(script, sep="\n", fill=TRUE)
      sink(con, type="message")
      e <- try(eval(parse(text=script), envir=env), silent=TRUE)
      sink(NULL, type="message")
      close(con)
      
      if (nzchar(str)) {
        info.bar$setMessageType(GtkMessageType["warning"])
        str <- paste(str, collapse="\n")
        str <- strsplit(str, " ")[[1]]
        str <- capture.output(cat(str, fill=80))
        str <- paste(str, collapse="\n")
        info.bar$setText(str)
        info.bar$show()
        p.current <<- e$plot
        p.list    <<- c(list(p.current), p.list)
        p.current.num <<- 1
        if(length(p.list) >= 2) button.prev["sensitive"] <<- TRUE
        button.next["sensitive"] <<- FALSE
      } else if (!is.list(e)) {
        info.bar$setMessageType(GtkMessageType["error"])
        str <- paste(e[1], collapse="\n")
        str <- strsplit(str, " ")[[1]]
        str <- capture.output(cat(str, fill=80))
        str <- paste(str, collapse="\n")
        info.bar$setText(str)
        info.bar$show()
      } else {
        p.current <<- e$plot
        p.list    <<- c(list(p.current), p.list)
        p.current.num <<- 1
        if(length(p.list) >= 2) button.prev["sensitive"] <<- TRUE
        button.next["sensitive"] <<- FALSE
        info.bar$hide()
      }
      
    },
        
    setModel = function(model){
      model <<- model
      if (constructed) {
        geom$completionSetModel(model)
        stratum$completionSetModel(model)
      }
    },
    
    setData = function(data){
      data <<- data
      rzPlotScript$setData(data)
      if (constructed) {
        misc$generateScript()
      }
    },
    
    onDetach = function(button){
      if(is.null(parent)) {
        main$hide()
        if (rzSettings$getEmbededDeviceOn()) {
          main$getChild1()$destroy()
          try(dev.off(), silent=TRUE)
        }
        button["label"] <- gettext("Attach")
        win <- gtkWindowNew(show=FALSE)
        parent.win <<- win
        if (rzSettings$getEmbededDeviceOn()) {
          win$setDefaultSize(900, 500)
        } else {
          win$setDefaultSize(400, 500)
        }
        win$setTitle(gettext("Plot"))
        parent <<- main$getParent()
        main$setOrientation(GtkOrientation["horizontal"])
        main$reparent(win)
        .self$construct()
        win$showAll()
        while (gtkEventsPending()) {}
        main$setPosition(500)
        while (gtkEventsPending()) {}
        if (rzSettings$getEmbededDeviceOn()) {
          p <- ggplot()
          p <- p + annotate("text", label = "Plot Area", x=10, y=10, size=14)
          p <- p + theme(axis.ticks=element_blank(), axis.text=element_blank(),
                         axis.title=element_blank())
          print(p)
        }
        
        gSignalConnect(win, "destroy", function(...){
          if (any(class(parent)=="GtkWidget")) {
            main$hide()
            if (rzSettings$getEmbededDeviceOn()) {
              main$getChild1()$destroy()
              try(dev.off(), silent=TRUE)
            }
            button["label"] <- gettext("Detach")
            main$setOrientation(GtkOrientation["vertical"])
            main$reparent(parent)
            parent.win <<- rzTools$getWindow()
            .self$construct()
            while (gtkEventsPending()) {}
            main$setPosition(300)
            while (gtkEventsPending()) {}
            if (rzSettings$getEmbededDeviceOn()) {
              p <- ggplot()
              p <- p + annotate("text", label = "Plot Area", x=10, y=10, size=14)
              p <- p + theme(axis.ticks=element_blank(), axis.text=element_blank(),
                             axis.title=element_blank())
              print(p)
            }
            parent <<- NULL            
          }
        })
        
      } else {
        main$getParent()$destroy()
      }
    }
))
rzplot$accessors("main", "constructed")
