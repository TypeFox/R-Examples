rzPlotLayerOptions <-
setRefClass("RzPlotLayerOptions",
  fields = c("layer", "type", "args", "mapping", "mapping.labels", "main","model"),
  method = list(
    initialize = function(...){
      initFields(...)
      args <<- list()
      mapping <<- list()
      mapping.labels <<- list()
      main <<- gtkDialogNewWithButtons(gettext("Options"), parent=rzTools$getWindow(),
                                       flags= c("modal", "destroy-with-parent"), 
                                       "gtk-close", GtkResponseType["accept"], show=FALSE)
      
      x     <- buildAes("x", combolabel=FALSE)
      y     <- buildAes("y", combolabel=FALSE)
      group <- buildAes("group", combolabel=FALSE)
      fill  <- buildAes("fill", combolabel=FALSE)
      color <- buildAes("colour", combolabel=FALSE)
      shape <- buildAes("shape", combolabel=FALSE)
      size  <- buildAes("size", combolabel=FALSE)
      line  <- buildAes("line", combolabel=FALSE)
      
      table.aes  <- gtkTableNew(homogeneous=FALSE)
      table.aes["border-width"] <- 5
      table.aes$setColSpacings(5)
      table.aes$setRowSpacings(2)
      table.aes$attach(x[[1]]    , 0, 1,  0,  1, "shrink", "shrink")
      table.aes$attach(x[[2]]    , 1, 2,  0,  1, 5       , "shrink")
#      table.aes$attach(x[[3]]    , 2, 3,  0,  1, 5       , "shrink")
      table.aes$attach(y[[1]]    , 0, 1,  1,  2, "shrink", "shrink")
      table.aes$attach(y[[2]]    , 1, 2,  1,  2, 5       , "shrink")
#      table.aes$attach(y[[3]]    , 2, 3,  1,  2, 5       , "shrink")
      table.aes$attach(group[[1]], 0, 1,  2,  3, "shrink", "shrink")
      table.aes$attach(group[[2]], 1, 2,  2,  3, 5       , "shrink")
#      table.aes$attach(group[[3]], 2, 3,  2,  3, 5       , "shrink")
      table.aes$attach(fill[[1]] , 0, 1,  3,  4, "shrink", "shrink")
      table.aes$attach(fill[[2]] , 1, 2,  3,  4, 5       , "shrink")
#      table.aes$attach(fill[[3]] , 2, 3,  3,  4, 5       , "shrink")
      table.aes$attach(color[[1]], 0, 1,  4,  5, "shrink", "shrink")
      table.aes$attach(color[[2]], 1, 2,  4,  5, 5       , "shrink")
#      table.aes$attach(color[[3]], 2, 3,  4,  5, 5       , "shrink")
      table.aes$attach(shape[[1]], 0, 1,  5,  6, "shrink", "shrink")
      table.aes$attach(shape[[2]], 1, 2,  5,  6, 5       , "shrink")
#      table.aes$attach(shape[[3]], 2, 3,  5,  6, 5       , "shrink")
      table.aes$attach(size[[1]] , 0, 1,  6,  7, "shrink", "shrink")
      table.aes$attach(size[[2]] , 1, 2,  6,  7, 5       , "shrink")
#      table.aes$attach(size[[3]] , 2, 3,  6,  7, 5       , "shrink")
      table.aes$attach(line[[1]] , 0, 1,  7,  8, "shrink", "shrink")
      table.aes$attach(line[[2]] , 1, 2,  7,  8, 5       , "shrink")
#      table.aes$attach(line[[3]] , 2, 3,  7,  8, 5       , "shrink")
      
      frame.aes <- gtkFrameNew("Layer Specific Variables")
      frame.aes$setShadowType(GtkShadowType["etched-in"])
      frame.aes$add(table.aes)
      
      # Appearance Options
      linetypes <- c("blank", "solid", "dashed", "dotted",
                     "dotdash", "longdash", "twodash")
      vbox.options <- gtkVBoxNew()
      fill2  <- .self$buildComboEntry("fill", colors(), width=120)
      color2 <- .self$buildComboEntry("color", colors(), width=120)
      shape2 <- .self$buildEntry("shape", class="any", width=120)
      size2  <- .self$buildEntry("size", class="numeric", width=120)
      width2 <- .self$buildEntry("width", class="numeric", width=120)
      line2  <- .self$buildCombo("linetype", linetypes, active=1, width=120)
      
      table.option1  <- gtkTableNew(homogeneous=FALSE)
      table.option1["border-width"] <- 5
      table.option1$setColSpacings(5)
      table.option1$setRowSpacings(2)
      table.option1$attach(fill2[[1]] , 0, 1, 0, 1, "shrink", "shrink")
      table.option1$attach(fill2[[2]] , 1, 2, 0, 1, 5       , "shrink")
      table.option1$attach(color2[[1]], 0, 1, 1, 2, "shrink", "shrink")
      table.option1$attach(color2[[2]], 1, 2, 1, 2, 5       , "shrink")
      table.option1$attach(shape2[[1]], 0, 1, 2, 3, "shrink", "shrink")
      table.option1$attach(shape2[[2]], 1, 2, 2, 3, 5       , "shrink")
      table.option1$attach(size2[[1]] , 0, 1, 3, 4, "shrink", "shrink")
      table.option1$attach(size2[[2]] , 1, 2, 3, 4, 5       , "shrink")
      table.option1$attach(width2[[1]], 0, 1, 4, 5, "shrink", "shrink")
      table.option1$attach(width2[[2]], 1, 2, 4, 5, 5       , "shrink")
      table.option1$attach(line2[[1]] , 0, 1, 5, 6, "shrink", "shrink")
      table.option1$attach(line2[[2]] , 1, 2, 5, 6, 5       , "shrink")
      
      frame.option1 <- gtkFrameNew("Appearance Options")
      frame.option1$setShadowType(GtkShadowType["etched-in"])
      frame.option1$add(table.option1)
      
      # position
      positions <- c("default", "identity", "dodge", "fill", "stack", "jitter")
      position <- .self$buildCombo(gettext("position"), positions, active=0, signal=FALSE)
      
      width.position  <-  .self$buildEntry("width", class="numeric", signal=FALSE)
      height.position <-  .self$buildEntry("height", class="numeric", signal=FALSE)
      
      table.position  <- gtkTableNew(homogeneous=FALSE)
      table.position["border-width"] <- 5
      table.position$setColSpacings(5)
      table.position$setRowSpacings(2)
      table.position$attach(position[[1]]       , 0, 1, 0, 1, "shrink", "shrink")
      table.position$attach(position[[2]]       , 1, 2, 0, 1, 5       , "shrink")
      table.position$attach(width.position[[1]] , 0, 1, 1, 2, "shrink", "shrink")
      table.position$attach(width.position[[2]] , 1, 2, 1, 2, 5       , "shrink")
      table.position$attach(height.position[[1]], 0, 1, 2, 3, "shrink", "shrink")
      table.position$attach(height.position[[2]], 1, 2, 2, 3, 5       , "shrink")
      frame.position <- gtkFrameNew("Position")
      frame.position$setShadowType(GtkShadowType["etched-in"])
      frame.position$add(table.position)
      
      
      setPosition <- function(...){
        position.type <- localize(position[[2]]$getActiveText())
        if (position.type=="default") {
          args[["position"]] <<- NULL
        } else {
          width  <- check.class(localize(width.position[[2]]$getText()) , "numeric")
          height <- check.class(localize(height.position[[2]]$getText()), "numeric")
          size <- paste(c(sprintf("width=%s", width), sprintf("height=%s", height)), collapse=",")
          args[["position"]] <<- sprintf("position_%s(%s)", position.type, size)
        }
      }
      gSignalConnect(position[[2]]       , "changed", setPosition)
      gSignalConnect(width.position[[2]] , "changed", setPosition)
      gSignalConnect(height.position[[2]], "changed", setPosition)
      
      # Options
      frame.option2 <- .self$buildOptions()
      
      # packing
      vbox <- main$getContentArea()
      hbox <- gtkHBoxNew(spacing=5)
      hbox$setBorderWidth(5)
      vbox$packStart(hbox, fill=TRUE)
      
      vbox2 <- gtkVBoxNew(spacing=5)
      vbox2$packStart(frame.option2, fill=TRUE)
      vbox2$packStart(frame.position, fill=TRUE)
      
      vbox3 <- gtkVBoxNew(spacing=5)
      vbox3$packStart(frame.option1 , fill=TRUE)
      
      vbox4 <- gtkVBoxNew(spacing=5)
      vbox4$packStart(frame.aes, fill=TRUE)

      hbox$packStart(vbox2, fill=TRUE)
      hbox$packStart(vbox3, fill=TRUE)
      hbox$packStart(vbox4, fill=TRUE)
      
      gSignalConnect(main, "response", function(dialog, respons.id){
        main$hideAll()
      })
    },
    
    getOptions = function(){
      main$showAll()
      main$run()
      if (length(mapping)==0) {
        return(args)        
      } else {
        aes <- c("x", "y", "group", "fill", "colour", "shape", "size", "line" )
        mapping.tmp <- unlist(mapping[aes])
        mapping.tmp <- sprintf("%s=%s", names(mapping.tmp), mapping.tmp)
        mapping.tmp <- paste(mapping.tmp, collapse=",")
        if (nzchar(mapping.tmp)) {
          mapping.tmp <- sprintf("aes(%s)", mapping.tmp)
          return(c(list(mapping.tmp), args))
        } else {
          return(args)
        }
        
      }
    },
    
    getMappingLabels = function(){
      aes <- c("x", "y", "group", "fill", "colour", "shape", "size", "linetype" )
      mapping.tmp <- unlist(mapping[aes])
      mapping.labels.tmp <- unlist(mapping.labels[names(mapping.tmp)])
      mapping.labels.tmp <- as.data.frame(rbind(vars=mapping.tmp, labels=mapping.labels.tmp), stringsAsFactors=FALSE)
      return(mapping.labels.tmp)
    },
    
    buildAes = function(label, combolabel=TRUE) {
      
      label.obj <- gtkLabelNew(label)
      combo <- new("RzCompletionCombo", width=150, show.combo=FALSE)
      
      gSignalConnect(combo$getCombo(), "changed", function(...){
        text <- localize(combo$getActiveText())
        if (nzchar(text)) {
          mapping[[label]] <<- text
        } else {
          mapping[[label]] <<- NULL          
        }
      })
      
      gSignalConnect(combo$getCombo(), "show", function(...){
        combo$setModel(model)
      })
      
      if (combolabel) {
        labels <- c(gettext("variable label"), gettext("variable name"), gettext("(free text)"))
        combo.label <- gtkComboBoxEntryNewText()
        for(i in labels) combo.label$appendText(i)
        entry.completion <- gtkEntryCompletionNew()
        entry.completion$setTextColumn(0)
        entry.completion$setInlineCompletion(TRUE)
        entry.completion$setInlineSelection(TRUE)
        entry.completion$setPopupSetWidth(FALSE)
        model.combo <- rGtkDataFrameNew(data.frame(labels))
        entry.completion$setModel(model.combo)
        entry <- combo.label$getChild()
        entry$setCompletion(entry.completion)
        
        gSignalConnect(combo.label, "changed", function(...){
          text <- localize(combo.label$getActiveText())
          mapping.labels[[label]] <<- text
        })
        
        combo.label$setActive(0)

        return(list(label.obj, combo$getCombo(), combo.label))
      } else {
        return(list(label.obj, combo$getCombo()))
      }
      
    },
    
    buildCombo = function(label=NULL, contents, active=NULL, class="character", width=120, signal=TRUE){
      combo <- gtkComboBoxNewText()
      combo$setData("class", class)
      for(i in contents) combo$appendText(i)
      combo["width-request"] <- width
      if (!is.null(active)) combo$setActive(active)
      default <- localize(combo$getActiveText())
      default <- check.class(default, class)
      
      if (signal) {
        gSignalConnect(combo, "changed", function(combo){
          text <- localize(combo$getActiveText())
          text <- check.class(text, class)
          if ( is.null(text) || (!is.null(default) && text==default) ) {
            args[[label]] <<- NULL
          } else {
            args[[label]] <<- deparse(text)
          }
        })
      }
      
      if (is.null(label)) {
        return(combo)
      } else {
        label.obj <- gtkLabelNew(label)
        return(list(label.obj, combo))        
      }
      
    },
    
    buildComboEntry = function(label=NULL, contents, active=NULL, class="character", width=120, signal=TRUE){
      combo <- gtkComboBoxEntryNewText()
      combo$setData("class", class)
      for(i in contents) combo$appendText(i)
      entry.completion <- gtkEntryCompletionNew()
      entry.completion$setTextColumn(0)
      entry.completion$setInlineCompletion(TRUE)
      entry.completion$setInlineSelection(TRUE)
      entry.completion$setPopupSetWidth(FALSE)
      model.combo <- rGtkDataFrameNew(data.frame(contents))
      entry.completion$setModel(model.combo)
      entry <- combo$getChild()
      entry$setCompletion(entry.completion)
      combo["width-request"] <- width
      if (!is.null(active)) combo$setActive(active)
      default <- localize(combo$getActiveText())
      default <- check.class(default, class)
      
      if (signal) {
        gSignalConnect(combo, "changed", function(combo){
          text <- localize(combo$getActiveText())
          text <- check.class(text, class)
          if ( is.null(text) || (!is.null(default) && text==default) ) {
            args[[label]] <<- NULL
          } else {
            args[[label]] <<- deparse(text)
          }
        })
      }
      
      if (is.null(label)) {
        return(combo)
      } else {
        label.obj <- gtkLabelNew(label)
        return(list(label.obj, combo))
      }
    },
    
    buildEntry = function(label, class="character", width=120, signal=TRUE){
      entry <- gtkEntryNew()
      entry["width-request"] <- width
      entry$setData("class", class)
      
      default <- localize(entry$getText())
      default <- check.class(default, class)
      
      if (signal) {
        gSignalConnect(entry, "changed", function(entry){
          text <- localize(entry$getText())
          text <- check.class(text, class)
          if ( is.null(text) || (!is.null(default) && text==default) ) {
            args[[label]] <<- NULL
          } else {
            args[[label]] <<- deparse(text)
          }
        })
      }
            
      label.obj <- gtkLabelNew(label)
      return(list(label.obj, entry))
    },
    
    buildOptions = function(){
      options <- list()
      if (layer=="geom") {
        if (type=="histogram" | type=="bar" | type=="freqpoly") {
          options[[1]] <- .self$buildEntry("width", class="numeric")
          options[[2]] <- .self$buildCombo("drop" , c("TRUE", "FALSE"), 1, class="logical")
          options[[3]] <- .self$buildCombo("right" , c("TRUE", "FALSE"), 1, class="logical")
          options[[4]] <- .self$buildEntry("binwidth", class="numeric")
          options[[5]] <- .self$buildEntry("origin", class="numeric")
          options[[6]] <- .self$buildEntry("break", class="caracter")
        } else if (type=="point") {
          options[[1]] <- .self$buildCombo("na.rm" , c("TRUE", "FALSE"), 1, class="logical")        
        } else if (type=="density") {
          options[[1]] <- .self$buildCombo("na.rm" , c("TRUE", "FALSE"), 1, class="logical")
          options[[2]] <- .self$buildEntry("adjust", class="numeric")
          kernel <-  c("gaussian", "epanechnikov", "rectangular",
                       "triangular", "biweight", "cosine", "optcosine")
          options[[3]] <- .self$buildCombo("kernel", kernel, 0, class="character")
          options[[4]] <- .self$buildCombo("trim" , c("TRUE", "FALSE"), 1, class="logical")
          
        } else if (type=="boxplot") {
          options[[1]] <- .self$buildComboEntry("outlier.colour", colors(), which(colors()=="black")-1, class="character")
          options[[2]] <- .self$buildEntry("outlier.shape", class="any")
          options[[3]] <- .self$buildEntry("outlier.size" , class="numeric")
          options[[4]] <- .self$buildCombo("notch" , c("TRUE", "FALSE"), 1, class="logical")
          options[[5]] <- .self$buildEntry("notchwidth" , class="numeric")
          options[[6]] <- .self$buildCombo("na.rm" , c("TRUE", "FALSE"), 1, class="logical")
          options[[7]] <- .self$buildEntry("coef" , class="numeric")
          
        } else if (type=="violin") {
          options[[1]] <- .self$buildCombo("trim" , c("TRUE", "FALSE"), 0, class="logical")
          options[[2]] <- .self$buildCombo("scale", c("equal", "count"), 0, class="logical")
          options[[3]] <- .self$buildEntry("adjust", class="numeric")
          kernel <-  c("gaussian", "epanechnikov", "rectangular",
                       "triangular", "biweight", "cosine", "optcosine")
          options[[4]] <- .self$buildCombo("kernel", kernel, 0, class="character")
          options[[5]] <- .self$buildCombo("na.rm" , c("TRUE", "FALSE"), 1, class="logical")
          
        } else if (type=="area" | type=="point" | type=="jitter") {
          options[[1]] <- .self$buildCombo("na.rm" , c("TRUE", "FALSE"), 1, class="logical")          
        } else if (type=="dotplot") {
          options <- list(
            .self$buildCombo("na.rm" , c("TRUE", "FALSE"), 1, class="logical"),
            .self$buildEntry("binwidth", class="numeric"),
            .self$buildCombo("binaxis" , c("x", "y"), 0, class="character"),
            .self$buildCombo("method" , c("dotdensity", "histodot"), 0, class="character"),
            .self$buildCombo("binpositions" , c("bygroup", "all"), 0, class="character"),
            .self$buildCombo("dir" , c("up", "down", "center", "centerwhole"), 0, class="character"),
            .self$buildEntry("dotsize", class="numeric"),
            .self$buildCombo("stackgroups" , c("TRUE", "FALSE"), 1, class="logical")
          )
        } else if (type=="rug") {
          options <- list(
            .self$buildEntry("sides" , class="character")
          )
        } else if (type=="smooth") {
          options <- list(
            .self$buildCombo("method" , c("auto", "lm", "glm", "gam", "loess", "rlm"), 0, class="character"),
            .self$buildEntry("formula" , class="formula"),            
            .self$buildCombo("se" , c("TRUE", "FALSE"), 0, class="logical"),
            .self$buildEntry("n" , class="numeric"),            
            .self$buildCombo("fullrange" , c("TRUE", "FALSE"), 1, class="logical"),
            .self$buildEntry("level" , class="numeric"),            
            .self$buildCombo("na.rm" , c("TRUE", "FALSE"), 1, class="logical")
          )
          
        } else if (type=="quantile") {
          options <- list(
            .self$buildCombo("lineend" , c("round", "butt", "square"), 1, class="character"),
            .self$buildCombo("linejoin", c("round", "mitre", "bevel"), 0, class="character"),
            .self$buildEntry("linemitre" , class="numeric"),            
            .self$buildCombo("na.rm" , c("TRUE", "FALSE"), 1, class="logical"),
            .self$buildEntry("quantiles", class="numeric"),            
            .self$buildEntry("formula" , class="formula"),            
            .self$buildCombo("method", c("rq"), 0, class="character")
          )
        }
      } else if (layer=="stat") {
        if (type=="summary") {
          fun <- c("none", "mean_se", "mean_sdl", "mean_cl_normal", "mean_cl_boot", "median_hilow")
          geoms <- c("errorbar", "pointrange", "linerange", "crossbar", "smooth", "point", "line", "bar")
          options <- list(
            .self$buildComboEntry("fun.data" , fun, 0, class="character"),
            .self$buildEntry("fun.y", class="character"),
            .self$buildEntry("fun.ymax", class="character"),
            .self$buildEntry("fun.ymin", class="character"),
            .self$buildComboEntry("geom", geoms, 1, class="character")
          )
        } else if (type=="qq") {
          options <- list(
            .self$buildEntry("distribution", class="character"),
            .self$buildEntry("dparams", class="list"),
            .self$buildCombo("na.rm" , c("TRUE", "FALSE"), 1, class="logical")
          )
        }
      }
      table <- gtkTableNew()
      table["border-width"] <- 5
      table$setColSpacings(5)
      table$setRowSpacings(2)
      for (i in seq_along(options)) {
        table$attach(options[[i]][[1]], 0, 1, i-1, i, "shrink", "shrink")
        table$attach(options[[i]][[2]], 1, 2, i-1, i, 5       , "shrink")
      }
      if (length(table$getChildren())==0) {
        table$attachDefaults(gtkLabelNew(gettext("No Option")), 0, 1, 0, 1)
      }
      frame <- gtkFrameNew("Layer Options")
      frame$setShadowType(GtkShadowType["etched-in"])
      frame$add(table)
      return(frame)
    }
  )
)
rzPlotLayerOptions$accessors(c("main", "args", "type", "model"))
