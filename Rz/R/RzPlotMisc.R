rzplot.misc <- 
setRefClass("RzPlotMisc",
  fields = c("axisPage", "themePage", "rzPlotScript",
             "flip.togglebutton",
             "combo.scalex", "combo.scaley",
             "combo.coordx", "combo.coordy",
             "entry.xlim.scale", "entry.ylim.scale",
             "entry.xlim.coord", "entry.ylim.coord",
             "title.entry", "xlab.combo", "ylab.combo"),
  methods = list(
    initialize  = function(...) {
      initFields(...)
      
      # misc options
      label.scale <- gtkLabelNew("trans")
      label.coord <- gtkLabelNew("trans")
      
      label.scalex <- gtkLabelNew("x")
      combo.scalex <<- gtkComboBoxNewText()
      scales <- c("default","log10", "reverse", "sqrt")
      for(i in scales) combo.scalex$appendText(i)
      combo.scalex$setActive(0)
      
      label.scaley <- gtkLabelNew("y")
      combo.scaley <<- gtkComboBoxNewText()
      for(i in scales) combo.scaley$appendText(i)
      combo.scaley$setActive(0)
      
      
      label.coordx <- gtkLabelNew("x")
      trans <- c("identity", "asn", "atanh","exp", "log",         
                 "log10", "log2", "logit", "probit", "reverse", "sqrt")
      combo.coordx <<- gtkComboBoxNewText()
      for(i in trans) combo.coordx$appendText(i)
      combo.coordx$setActive(0)
      
      label.coordy <- gtkLabelNew("y")
      combo.coordy <<- gtkComboBoxNewText()
      for(i in trans) combo.coordy$appendText(i)
      combo.coordy$setActive(0)
      
      label.limits.scale <-  gtkLabelNew(gettext("range"))
      entry.xlim.scale   <<- gtkEntryNew()
      entry.ylim.scale   <<- gtkEntryNew()
      entry.xlim.scale["width-request"] <<- 1
      entry.ylim.scale["width-request"] <<- 1

      label.limits.coord <-  gtkLabelNew(gettext("range"))
      entry.xlim.coord   <<- gtkEntryNew()
      entry.ylim.coord   <<- gtkEntryNew()
      entry.xlim.coord["width-request"] <<- 1
      entry.ylim.coord["width-request"] <<- 1
      
      #flip
      flip.togglebutton <<- gtkToggleButtonNewWithLabel(gettext("flip"))
      
      table.scale <- gtkTableNew(FALSE)
      table.scale["border-width"] <- 5
      table.scale$attach        (label.scalex      , 1, 2, 0, 1, "fill"  , "shrink", 0, 0)
      table.scale$attach        (label.scaley      , 2, 3, 0, 1, "fill"  , "shrink", 0, 0)
      table.scale$attach        (label.scale       , 0, 1, 1, 2, "shrink", "shrink", 0, 0)
      table.scale$attachDefaults(combo.scalex      , 1, 2, 1, 2)
      table.scale$attachDefaults(combo.scaley      , 2, 3, 1, 2)
      table.scale$attach        (label.limits.scale, 0, 1, 2, 3, "shrink", "shrink", 0, 0)
      table.scale$attachDefaults(entry.xlim.scale  , 1, 2, 2, 3)
      table.scale$attachDefaults(entry.ylim.scale  , 2, 3, 2, 3)
      table.scale$setColSpacings(5)
      table.scale$setRowSpacings(2)
      frame.scale <- gtkFrameNew("Scale")
      frame.scale$setShadowType(GtkShadowType["etched-in"])
      frame.scale$add(table.scale)
      
      table.coord <- gtkTableNew(FALSE)
      table.coord["border-width"] <- 5
      table.coord$attach        (label.coordx      , 1, 2, 0, 1, "fill"  , "shrink", 0, 0)
      table.coord$attach        (label.coordy      , 2, 3, 0, 1, "fill"  , "shrink", 0, 0)
      table.coord$attach        (label.coord       , 0, 1, 1, 2, "shrink", "shrink", 0, 0)
      table.coord$attachDefaults(combo.coordx      , 1, 2, 1, 2)
      table.coord$attachDefaults(combo.coordy      , 2, 3, 1, 2)
      table.coord$attach        (label.limits.coord, 0, 1, 2, 3, "shrink", "shrink", 0, 0)
      table.coord$attachDefaults(entry.xlim.coord  , 1, 2, 2, 3)
      table.coord$attachDefaults(entry.ylim.coord  , 2, 3, 2, 3)
      table.coord$attachDefaults(flip.togglebutton , 0, 3, 3, 4)
      table.coord$setColSpacings(5)
      table.coord$setRowSpacings(2)
      frame.coord <- gtkFrameNew("Coordinate")
      frame.coord$setShadowType(GtkShadowType["etched-in"])
      frame.coord$add(table.coord)
      
      vbox.axisPage <- gtkVBoxNew()
      vbox.axisPage$packStart(frame.scale, expand=FALSE)
      vbox.axisPage$packStart(frame.coord, expand=FALSE)
      axisPage <<- buildPlotOptionPage(vbox.axisPage)
      
      #theme
      button.customize <- gtkButtonNewWithLabel(gettext("Open Theme Editor"))
      hbox.button.customize <- gtkHBoxNew()
      hbox.button.customize$packStart(button.customize, expand=TRUE, fill=FALSE)
      gSignalConnect(button.customize, "clicked", function(...) {
        if (is.null(rzTools$getThemeEditor())) {
          rzTools$setThemeEditor(new("RzPlotTheme"))        
        }
        rzTools$getThemeEditor()$show()
      })
      
      table.theme <- gtkTableNew(FALSE)
      table.theme["border-width"] <- 5
      table.theme$setColSpacings(5)
      table.theme$setRowSpacings(2)
      table.theme$attachDefaults(hbox.button.customize, 0, 1, 0, 1)
      
      frame.theme <- gtkFrameNew("Theme")
      frame.theme$setShadowType(GtkShadowType["etched-in"])
      frame.theme$add(table.theme)
      
      vbox.themePage <- gtkVBoxNew()
      vbox.themePage$packStart(frame.theme, expand=FALSE)
      themePage <<- buildPlotOptionPage(vbox.themePage)

      .self$generateScript()
      
      gSignalConnect(combo.scalex      , "changed", .self$generateScript)
      gSignalConnect(combo.scaley      , "changed", .self$generateScript)
      gSignalConnect(combo.coordx      , "changed", .self$generateScript)
      gSignalConnect(combo.coordy      , "changed", .self$generateScript)
      gSignalConnect(entry.xlim.scale  , "changed", .self$generateScript)
      gSignalConnect(entry.ylim.scale  , "changed", .self$generateScript)
      gSignalConnect(entry.xlim.coord  , "changed", .self$generateScript)
      gSignalConnect(entry.ylim.coord  , "changed", .self$generateScript)
      gSignalConnect(flip.togglebutton , "toggled", .self$generateScript)
    },
    
    clear = function(){
      combo.scalex$setActive(0)
      combo.scaley$setActive(0)
      combo.coordx$setActive(0)
      combo.coordy$setActive(0)
      entry.xlim.scale$setText("")
      entry.ylim.scale$setText("")
      entry.xlim.coord$setText("")
      entry.ylim.coord$setText("")
      flip.togglebutton$setActive(FALSE)
    },
    
    generateScript = function(...){
      scalex <- localize(combo.scalex$getActiveText())
      scaley <- localize(combo.scaley$getActiveText())
      coordx <- localize(combo.coordx$getActiveText())
      coordy <- localize(combo.coordy$getActiveText())
      xlim.scale   <- localize(entry.xlim.scale$getText())
      ylim.scale   <- localize(entry.ylim.scale$getText())
      xlim.scale   <- strsplit(xlim.scale, ",")[[1]]
      ylim.scale   <- strsplit(ylim.scale, ",")[[1]]
      xlim.scale   <- suppressWarnings(as.numeric(xlim.scale))
      ylim.scale   <- suppressWarnings(as.numeric(ylim.scale))
      xlim.coord   <- localize(entry.xlim.coord$getText())
      ylim.coord   <- localize(entry.ylim.coord$getText())
      xlim.coord   <- strsplit(xlim.coord, ",")[[1]]
      ylim.coord   <- strsplit(ylim.coord, ",")[[1]]
      xlim.coord   <- suppressWarnings(as.numeric(xlim.coord))
      ylim.coord   <- suppressWarnings(as.numeric(ylim.coord))
      
      if (coordx == "identity") coordx <- NULL
      if (coordy == "identity") coordy <- NULL
      
      if (any(is.na(xlim.scale)) | length(xlim.scale) != 2) xlim.scale <- NULL
      if (any(is.na(ylim.scale)) | length(ylim.scale) != 2) ylim.scale <- NULL
      if (any(is.na(xlim.coord)) | length(xlim.coord) != 2) xlim.coord <- NULL
      if (any(is.na(ylim.coord)) | length(ylim.coord) != 2) ylim.coord <- NULL
      
      flip   <- flip.togglebutton$getActive()
      
      if (scalex != "default") rzPlotScript$setScript(layer="scale_x", type=scalex)
      else                     rzPlotScript$clearScript("scale_x")

      if (scaley != "default") rzPlotScript$setScript("scale_y", type=scaley)
      else                     rzPlotScript$clearScript("scale_y")
      
      if (!is.null(coordx) || !is.null(coordy)) {
        rzPlotScript$setScript("coord", type="trans",
                               args=list(xtrans=deparse(coordx), ytrans=deparse(coordy), limx=deparse(xlim), limy=deparse(ylim)))
        if (flip) rzPlotScript$setScript("coord", type="flip", add=TRUE)
        
      } else if (!is.null(xlim.coord) || !is.null(ylim.coord)) {
        rzPlotScript$setScript("coord", type="cartesian", args=list(xlim=deparse(xlim.coord), ylim=deparse(ylim.coord)))
        if (flip) rzPlotScript$setScript("coord", type="flip", add=TRUE)
        
      } else if (flip) {
        rzPlotScript$setScript("coord", type="flip")
        
      } else {
        rzPlotScript$clearScript("coord")
      }
      
      if (is.null(xlim.scale)) {
        rzPlotScript$setFreeScript(name="xlim", script=NA)
      } else {
        rzPlotScript$setFreeScript(name="xlim", script=sprintf("xlim(%s)", deparse(xlim.scale)))
      }
      
      if (is.null(ylim.scale)) {
        rzPlotScript$setFreeScript(name="ylim", script=NA)
      } else {
        rzPlotScript$setFreeScript(name="ylim", script=sprintf("ylim(%s)", deparse(ylim.scale)))
      }
      
    }
  )
)
rzplot.misc$accessors(c("axisPage", "themePage"))
