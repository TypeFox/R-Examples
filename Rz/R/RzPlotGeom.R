rzplot.geom <- 
setRefClass("RzPlotGeom",
  fields = c("main", "rzPlotScript",
             "combo.geom", "combo.stat", "options", "options.obj", "labels",
             "combo.x", "combo.y", "xlab.combo", "ylab.combo", "title.entry", "model"),
  methods = list(
    initialize  = function(...) {
      initFields(...)
      model <<- NULL
      options <<- list()
      options.obj <<- list()
      labels <<- list()
      labels[["geom"]] <<- data.frame(1:2, row.names=c("vars", "labels"), stringsAsFactors=FALSE)
      labels[["geom"]][[1]] <<- NULL
      labels[["stat"]] <<- data.frame(1:2, row.names=c("vars", "labels"), stringsAsFactors=FALSE)
      labels[["stat"]][[1]] <<- NULL
      
      # variables
      label.x <- gtkLabelNew(gettext("x"))
      combo.x <<- new("RzCompletionCombo", width=40)      
      xlab.combo <<- gtkComboBoxEntryNewText()
      xlab.combo["width-request"] <<- 1
      xlab.combo$show()
      labs <- c(gettext("variable label"), gettext("variable name"), gettext("(free text)"))
      for(i in labs) xlab.combo$appendText(i)
      xlab.combo$setActive(0)
      
      label.y <- gtkLabelNew(gettext("y"))
      combo.y <<- new("RzCompletionCombo", width=40)
      ylab.combo <<- gtkComboBoxEntryNewText()
      ylab.combo["width-request"] <<- 1
      ylab.combo$show()
      for(i in labs) ylab.combo$appendText(i)
      ylab.combo$setActive(0)
      
      title.label <-  gtkLabelNew(gettext("title"))
      title.entry <<- gtkEntryNew()
      
      frame.var <- gtkFrameNew()
      frame.var$setShadowType(GtkShadowType["etched-in"])
      table.var <- gtkTableNew(homogeneous=FALSE)
      table.var["border-width"] <- 5
      table.var$attach(label.x           , 0, 1, 0, 1, "shrink", "shrink", 5)
      table.var$attach(combo.x$getCombo(), 1, 2, 0, 1, 5       , "shrink")
      table.var$attach(xlab.combo        , 2, 3, 0, 1, 5       , "shrink")
      table.var$attach(label.y           , 0, 1, 1, 2, "shrink", "shrink", 5)
      table.var$attach(combo.y$getCombo(), 1, 2, 1, 2, 5       , "shrink")
      table.var$attach(ylab.combo        , 2, 3, 1, 2, 5       , "shrink")
      table.var$attach(title.label       , 0, 1, 2, 3, "shrink", "shrink", 5)
      table.var$attach(title.entry       , 1, 3, 2, 3, 5       , "shrink")
      table.var$setColSpacings(5)
      table.var$setRowSpacings(2)
      frame.var$add(table.var)
      
      # geom
      label.geom <- gtkLabelNew("geom")
      combo.geom <<- gtkComboBoxNewText()
      combo.geom["width-request"] <<- 1
      geoms <- c("none", "bar", "freqpoly", "histogram", "density", "boxplot", "violin",
                 "line", "area", "point", "jitter", "dotplot", "rug", "smooth", "quantile", "blank")
      for(i in geoms) combo.geom$appendText(i)
      combo.geom$setActive(0)
      
      button.geom.options <- gtkButtonNewWithLabel(gettext("Options"))

      image  <- gtkImageNewFromStock(GTK_STOCK_ADD, GtkIconSize["menu"])
      button.geom.stack <- gtkButtonNew()
      button.geom.stack["tooltip-text"] <- gettext("Stack Layer")
      button.geom.stack$setFocusOnClick(FALSE)
      button.geom.stack$setImage(image)
      button.geom.stack$setRelief(GtkReliefStyle["none"])
      
      image   <- gtkImageNewFromStock(GTK_STOCK_REMOVE, GtkIconSize["menu"])
      button.geom.remove <- gtkButtonNew()
      button.geom.remove["tooltip-text"] <- gettext("Remove Layer")
      button.geom.remove$setFocusOnClick(FALSE)
      button.geom.remove$setImage(image)
      button.geom.remove$setRelief(GtkReliefStyle["none"])
            
      # stat
      label.stat <- gtkLabelNew("stat")
      combo.stat <<- gtkComboBoxNewText()
      combo.stat["width-request"] <<- 1
      stats <- c("none", "sum", "summary", "qq")
      for(i in stats) combo.stat$appendText(i)
      combo.stat$setActive(0)
      
      button.stat.options <- gtkButtonNewWithLabel(gettext("Options"))
      
      image  <- gtkImageNewFromStock(GTK_STOCK_ADD, GtkIconSize["menu"])
      button.stat.stack <- gtkButtonNew()
      button.stat.stack["tooltip-text"] <- gettext("Stack Layer")
      button.stat.stack$setFocusOnClick(FALSE)
      button.stat.stack$setImage(image)
      button.stat.stack$setRelief(GtkReliefStyle["none"])
      
      image   <- gtkImageNewFromStock(GTK_STOCK_REMOVE, GtkIconSize["menu"])
      button.stat.remove <- gtkButtonNew()
      button.stat.remove["tooltip-text"] <- gettext("Remove Layer")
      button.stat.remove$setFocusOnClick(FALSE)
      button.stat.remove$setImage(image)
      button.stat.remove$setRelief(GtkReliefStyle["none"])
      
      # packing
      table.layer <- gtkTableNew(homogeneous=FALSE)
      table.layer["border-width"] <- 5
      table.layer$setColSpacings(5)
      table.layer$setRowSpacings(0)
      table.layer$attach        (label.geom         , 0, 1, 0, 1, "shrink", "shrink", 2)
      table.layer$attachDefaults(combo.geom         , 1, 2, 0, 1)
      table.layer$attach        (button.geom.options, 2, 3, 0, 1, "shrink", "shrink", 2)
      table.layer$attach        (button.geom.stack  , 3, 4, 0, 1, "shrink", "shrink")
      table.layer$attach        (button.geom.remove , 4, 5, 0, 1, "shrink", "shrink")
      table.layer$attach        (label.stat         , 0, 1, 1, 2, "shrink", "shrink", 2)
      table.layer$attachDefaults(combo.stat         , 1, 2, 1, 2)
      table.layer$attach        (button.stat.options, 2, 3, 1, 2, "shrink", "shrink", 2)
      table.layer$attach        (button.stat.stack  , 3, 4, 1, 2, "shrink", "shrink")
      table.layer$attach        (button.stat.remove , 4, 5, 1, 2, "shrink", "shrink")
      
      frame.layer <- gtkFrameNew("Layer")
      frame.layer$setShadowType(GtkShadowType["etched-in"])
      frame.layer$add(table.layer)
      
      main.vbox <- gtkVBoxNew()
      main.vbox$packStart(frame.var  , expand=FALSE)
      main.vbox$packStart(frame.layer, expand=FALSE)
      
      main <<- buildPlotOptionPage(main.vbox)

      gSignalConnect(combo.geom, "changed", function(...){
        options[["geom"]] <<- NULL
        options.obj[["geom"]] <<- NULL
        labels[["geom"]] <<- data.frame(1:2, row.names=c("vars", "labels"))
        labels[["geom"]][[1]] <<- NULL
        
        .self$generateScript()
      })
      gSignalConnect(combo.stat, "changed", function(...){
        options[["stat"]] <<- NULL
        options.obj[["stat"]] <<- NULL
        labels[["stat"]] <<- data.frame(1:2, row.names=c("vars", "labels"))
        labels[["stat"]][[1]] <<- NULL
        .self$generateScript()
      })
      
      gSignalConnect(combo.x$getCombo(), "changed", .self$generateScript)
      gSignalConnect(combo.y$getCombo(), "changed", .self$generateScript)
      gSignalConnect(xlab.combo        , "changed", .self$generateScript)
      gSignalConnect(ylab.combo        , "changed", .self$generateScript)
      gSignalConnect(title.entry       , "changed", .self$generateScript)
      
      gSignalConnect(button.geom.stack  , "clicked", .self$stackLayer , "geom")
      gSignalConnect(button.stat.stack  , "clicked", .self$stackLayer , "stat")
      gSignalConnect(button.geom.remove , "clicked", .self$removeLayer, "geom")
      gSignalConnect(button.stat.remove , "clicked", .self$removeLayer, "stat")
      gSignalConnect(button.geom.options, "clicked", .self$setOptions , "geom")
      gSignalConnect(button.stat.options, "clicked", .self$setOptions , "stat")
      
      .self$generateScript()
      
    },
    
    stackLayer = function(button, data){
      layer <- data[[1]]
      rzPlotScript$stackLayer(layer)
    },

    removeLayer = function(button, data){
      layer <- data[[1]]
      rzPlotScript$removeLayer(layer)
    },
    
    setOptions = function(button, data){
      layer <- data[[1]]
      type  <- NULL
      if (layer=="geom") {
        type <- combo.geom$getActiveText()
      } else {
        type <- combo.stat$getActiveText()        
      }
      
      if (type=="none") return()
      
      if (is.null(options.obj[[layer]]) || options.obj[[layer]]$getType() != type) {
        options.obj[[layer]] <<- new("RzPlotLayerOptions", layer=layer, type=type, model=model)
        
      } else {
        options.obj[[layer]]$setModel(model)
      }
      options[[layer]] <<- options.obj[[layer]]$getOptions()
#      labels[[layer]]  <<-  options.obj[[layer]]$getMappingLabels()
      
      .self$generateScript()
    },
    
    generateScript = function(...){
      geom  <- combo.geom$getActiveText()
      stat  <- combo.stat$getActiveText()
      x     <- localize(combo.x$getActiveText())
      y     <- localize(combo.y$getActiveText())
      xlab  <- localize(xlab.combo$getActiveText())
      ylab  <- localize(ylab.combo$getActiveText())
      title <- localize(title.entry$getText())
      
      if (geom=="none") {
        rzPlotScript$clearScript("geom")
      } else {
        rzPlotScript$setScript("geom", type=geom, args=options[["geom"]])
      }
      
      if (stat=="none") {
        rzPlotScript$clearScript("stat")
      } else {
        rzPlotScript$setScript("stat", type=stat, args=options[["stat"]])
      }
      
      rzPlotScript$setAes(c("x", "y"), c(x, y))
      
#      labels.tmp <- labels[["geom"]]
#      labels.tmp[colnames(labels[["stat"]])] <- labels[["stat"]]        

      labels.tmp        <- data.frame(1:2, row.names=c("vars", "labels"), stringsAsFactors=FALSE)
      labels.tmp[[1]]   <- NULL
      labels.tmp[["x"]] <- c(x , xlab)
      labels.tmp[["y"]] <- c(y , ylab)
      labels.tmp[["title"]] <- c("", title)
      rzPlotScript$setLabs(labels.tmp)
    },
    
    clear = function(){
      combo.geom$setActive(0)
      combo.stat$setActive(0)
      combo.x$clear()
      combo.y$clear()
      title.entry$setText("")
      xlab.combo$setActive(0)
      ylab.combo$setActive(0)
    },
    
    setX = function(txt){
      combo.x$setText(txt)
      .self$generateScript()
    },
    
    completionSetModel = function(model){
      model <<- model
      combo.x$setModel(model)
      combo.y$setModel(model)
      if (!is.null(options.obj[["geom"]])) options.obj[["geom"]]$setModel(model)
      if (!is.null(options.obj[["stat"]])) options.obj[["stat"]]$setModel(model)
    }
  )
)
rzplot.geom$accessors("main")
