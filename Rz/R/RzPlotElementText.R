rzplot.element.text <- 
setRefClass("RzPlotElementText",
  fields = c("main", "name", "family", "face", "colour", "size", "size.rel",
             "hjust", "vjust", "angle", "lineheight", "parent", "script", "blank",
             "family.button", "face.combo", "color.widget" ,"size.button", "size.rel.button",
             "hjust.button", "vjust.button", "angle.button", "lineheight.button",
             "inherit.from", "family.inherit", "face.inherit", "colour.inherit", "size.inherit",
             "hjust.inherit", "vjust.inherit", "angle.inherit", "lineheight.inherit",
             "inherit.family.button", "inherit.face.button", "inherit.colour.button",
             "inherit.size.button", "inherit.hjust.button", "inherit.vjust.button",
             "inherit.angle.button", "inherit.lineheight.button", "button.blank", "faces"),
  methods = list(
    initialize  = function(...) {
      initFields(...)
      family.inherit     <<- FALSE
      face.inherit       <<- FALSE
      colour.inherit     <<- FALSE
      size.inherit       <<- FALSE
      hjust.inherit      <<- FALSE
      vjust.inherit      <<- FALSE
      angle.inherit      <<- FALSE
      lineheight.inherit <<- FALSE
      if (class(family) == "uninitializedField"){
        family <<- "Sans"
        family.inherit <<- TRUE
      }
      if (class(face) == "uninitializedField"){
        face <<- "plain"
        face.inherit <<- TRUE
      }
      if (class(colour) == "uninitializedField") {
        colour <<- "#000000"
        colour.inherit <<- TRUE
      }
      if (class(size) == "uninitializedField") {
        size <<- 1
        size.rel <<- FALSE
        size.inherit <<- TRUE
      }
      if (class(hjust) == "uninitializedField"){
        hjust <<- 0
        hjust.inherit <<- TRUE
      }
      if (class(vjust) == "uninitializedField"){
        vjust <<- 0
        vjust.inherit <<- TRUE
      }
      if (class(angle) == "uninitializedField"){
        angle <<- 0
        angle.inherit <<- TRUE
      }
      if (class(lineheight) == "uninitializedField"){
        lineheight <<- 0
        lineheight.inherit <<- TRUE
      }
      if (class(blank) == "uninitializedField") {
        blank <<- FALSE
      }
      
      family.button <<- gtkFontSelectionButtonNew(fontname=family, parent=parent)
      faces <<-  c("plain", "italic", "bold", "bold.italic")
      face.combo <<- .self$buildCombo(faces, which(faces == face) - 1)
      color.widget <<- gtkColorSelectionWidgetNew(spacing=2, parent=parent)
      color.widget$setColor(colour)
      size.adj <- gtkAdjustmentNew(size, 0, 99, 0.1)
      size.button   <<- gtkSpinButtonNew(size.adj, climb.rate=0.1, digits=1)
      size.button$setValue(size)
      size.rel.button <<- gtkCheckButtonNewWithLabel("rel")
      size.rel.button$setActive(size.rel)
      size.rel.button["tooltip-text"] <<- gettext("Relative sizing")
      size.hbox <- gtkHBoxNew(spacing=2)
      size.hbox$packStart(size.button)
      size.hbox$packStart(size.rel.button, expand=FALSE)
      hjust.adj <- gtkAdjustmentNew(size, 0, 1, 0.01)
      hjust.button   <<- gtkSpinButtonNew(hjust.adj, climb.rate=0.01, digits=2)
      hjust.button$setValue(hjust)
      vjust.adj <- gtkAdjustmentNew(size, 0, 1, 0.01)
      vjust.button   <<- gtkSpinButtonNew(vjust.adj, climb.rate=0.01, digits=2)
      vjust.button$setValue(vjust)
      angle.adj <- gtkAdjustmentNew(size, -360, 360, 1)
      angle.button   <<- gtkSpinButtonNew(angle.adj, climb.rate=1, digits=0)
      angle.button$setValue(angle)
      lineheight.adj <- gtkAdjustmentNew(size, 0, 10, 0.1)
      lineheight.button   <<- gtkSpinButtonNew(lineheight.adj, climb.rate=0.1, digits=2)
      lineheight.button$setValue(lineheight)
      
      table <- gtkTableNew()
      table$setBorderWidth(5)
      table$attach(gtkLabelNew(gettext("Font\nFamily")), 0, 1, 0, 1, "shrink", "shrink")
      table$attach(family.button                      , 1, 2, 0, 1, 5  , "shrink")
      table$attach(gtkLabelNew(gettext("Font\nFace"))  , 0, 1, 1, 2, "shrink", "shrink")
      table$attach(face.combo                         , 1, 2, 1, 2, 5  , "shrink")
      table$attach(gtkLabelNew(gettext("Colour"))     , 0, 1, 2, 3, "shrink", "shrink")
      table$attach(color.widget                       , 1, 2, 2, 3, 5  , "shrink")
      table$attach(gtkLabelNew(gettext("Size"))       , 0, 1, 3, 4, "shrink", "shrink")
      table$attach(size.hbox                          , 1, 2, 3, 4, 5  , "shrink")
      table$attach(gtkLabelNew(gettext("Horizontal\nJustification")), 0, 1, 4, 5, "shrink", "shrink")
      table$attach(hjust.button                       , 1, 2, 4, 5, 5  , "shrink")
      table$attach(gtkLabelNew(gettext("Vertical\nJustification")), 0, 1, 5, 6, "shrink", "shrink")
      table$attach(vjust.button                       , 1, 2, 5, 6, 5  , "shrink")
      table$attach(gtkLabelNew(gettext("Angle"))      , 0, 1, 6, 7, "shrink", "shrink")
      table$attach(angle.button                       , 1, 2, 6, 7, 5  , "shrink")
      table$attach(gtkLabelNew(gettext("Line\nHeight")), 0, 1, 7, 8, "shrink", "shrink")
      table$attach(lineheight.button                  , 1, 2, 7, 8, 5  , "shrink")
      table$setRowSpacings(5)
      table$setColSpacings(2)
      
      button.blank <<- gtkToggleButtonNewWithLabel(name)
      button.blank["tooltip-markup"] <<- gettext("<span font_style='italic' font_weight='bold'>theme_blank()</span> when turn <span font_style='italic' font_weight='bold'>OFF</span> the button")
      button.blank$setActive(!blank)
      gSignalConnect(button.blank, "toggled", function(...){
        blank <<- ! button.blank$getActive()
        if(blank) {
          table$setSensitive(FALSE)          
        } else {
          table$setSensitive(TRUE)
        }
        .self$generateScript()
      })
      hbox.blank <- gtkHBoxNew(spacing=2)
      hbox.blank$packStart(button.blank, expand=FALSE)
      
      main <<- gtkFrameNew()
      main$setLabelWidget(hbox.blank)
      main$setShadowType(GtkShadowType["etched-in"])
      main$add(table)
      
      gSignalConnect(family.button, "clicked", function(obj){
        family   <<- family.button$getLabel()
        .self$generateScript()
      })
      gSignalConnect(face.combo, "changed", function(obj){
        face     <<- face.combo$getActiveText()
        .self$generateScript()
      })
      gSignalConnect(color.widget$getEntry(), "changed", function(obj){
        colour   <<- obj$getText()
        .self$generateScript()
      })
      gSignalConnect(size.button, "changed", function(obj){
        size     <<- suppressWarnings(as.numeric(obj$getText()))
        .self$generateScript()
      })
      gSignalConnect(size.rel.button, "toggled", function(obj){
        size.rel <<- obj$getActive()
        .self$generateScript()
      })
      gSignalConnect(hjust.button, "changed", function(obj){
        hjust    <<- suppressWarnings(as.numeric(obj$getText()))
        .self$generateScript()
      })
      gSignalConnect(vjust.button, "changed", function(obj){
        vjust    <<- suppressWarnings(as.numeric(obj$getText()))
        .self$generateScript()
      })
      gSignalConnect(angle.button, "changed", function(obj){
        angle    <<- suppressWarnings(as.numeric(obj$getText()))
        .self$generateScript()
      })
      gSignalConnect(lineheight.button, "changed", function(obj){
        lineheight <<- suppressWarnings(as.numeric(obj$getText()))
        .self$generateScript()
      })
      
      tooltip.text <- character()
      if (is.null(inherit.from)) {
        tooltip.text <- gettext("Inherit from the current theme")
      } else {
        tooltip.text <- gettextf("Inherit from <span font_style='italic' font_weight='bold' font_size='large'>%s</span>", inherit.from)
      }
      
      inherit.family.button     <<- gtkToggleButtonNew()
      inherit.face.button       <<- gtkToggleButtonNew()
      inherit.colour.button     <<- gtkToggleButtonNew()
      inherit.size.button       <<- gtkToggleButtonNew()
      inherit.hjust.button      <<- gtkToggleButtonNew()
      inherit.vjust.button      <<- gtkToggleButtonNew()
      inherit.angle.button      <<- gtkToggleButtonNew()
      inherit.lineheight.button <<- gtkToggleButtonNew()        
      inherit.family.button$setImage(gtkImageNewFromStock(GTK_STOCK_GO_DOWN, GtkIconSize["button"]))
      inherit.face.button$setImage(gtkImageNewFromStock(GTK_STOCK_GO_DOWN, GtkIconSize["button"]))
      inherit.colour.button$setImage(gtkImageNewFromStock(GTK_STOCK_GO_DOWN, GtkIconSize["button"]))
      inherit.size.button$setImage(gtkImageNewFromStock(GTK_STOCK_GO_DOWN, GtkIconSize["button"]))
      inherit.hjust.button$setImage(gtkImageNewFromStock(GTK_STOCK_GO_DOWN, GtkIconSize["button"]))
      inherit.vjust.button$setImage(gtkImageNewFromStock(GTK_STOCK_GO_DOWN, GtkIconSize["button"]))
      inherit.angle.button$setImage(gtkImageNewFromStock(GTK_STOCK_GO_DOWN, GtkIconSize["button"]))
      inherit.lineheight.button$setImage(gtkImageNewFromStock(GTK_STOCK_GO_DOWN, GtkIconSize["button"]))
      inherit.family.button["tooltip-markup"]     <<- tooltip.text
      inherit.face.button["tooltip-markup"]       <<- tooltip.text
      inherit.colour.button["tooltip-markup"]     <<- tooltip.text
      inherit.size.button["tooltip-markup"]       <<- tooltip.text
      inherit.hjust.button["tooltip-markup"]      <<- tooltip.text
      inherit.vjust.button["tooltip-markup"]      <<- tooltip.text
      inherit.angle.button["tooltip-markup"]      <<- tooltip.text
      inherit.lineheight.button["tooltip-markup"] <<- tooltip.text
      
      table$attach(inherit.family.button,     2, 3, 0, 1, "shrink", "shrink")
      table$attach(inherit.face.button,       2, 3, 1, 2, "shrink", "shrink")
      table$attach(inherit.colour.button,     2, 3, 2, 3, "shrink", "shrink")
      table$attach(inherit.size.button,       2, 3, 3, 4, "shrink", "shrink")
      table$attach(inherit.hjust.button,      2, 3, 4, 5, "shrink", "shrink")
      table$attach(inherit.vjust.button,      2, 3, 5, 6, "shrink", "shrink")
      table$attach(inherit.angle.button,      2, 3, 6, 7, "shrink", "shrink")
      table$attach(inherit.lineheight.button, 2, 3, 7, 8, "shrink", "shrink")
      
      gSignalConnect(inherit.family.button, "toggled", function(obj){
        family.inherit <<- obj$getActive()
        family.button$setSensitive( !family.inherit )
        .self$generateScript()
      })
      gSignalConnect(inherit.face.button, "toggled", function(obj){
        face.inherit <<- obj$getActive()
        face.combo$setSensitive( !face.inherit )
        .self$generateScript()
      })
      gSignalConnect(inherit.colour.button, "toggled", function(obj){
        colour.inherit <<- obj$getActive()
        color.widget$setSensitive( !colour.inherit )
        .self$generateScript()
      })
      gSignalConnect(inherit.size.button, "toggled", function(obj){
        size.inherit <<- obj$getActive()
        size.hbox$setSensitive( !size.inherit )
        .self$generateScript()
      })
      gSignalConnect(inherit.hjust.button, "toggled", function(obj){
        hjust.inherit <<- obj$getActive()
        hjust.button$setSensitive( !hjust.inherit )
        .self$generateScript()
      })
      gSignalConnect(inherit.vjust.button, "toggled", function(obj){
        vjust.inherit <<- obj$getActive()
        vjust.button$setSensitive( !vjust.inherit )
        .self$generateScript()
      })
      gSignalConnect(inherit.angle.button, "toggled", function(obj){
        angle.inherit <<- obj$getActive()
        angle.button$setSensitive( !angle.inherit )
        .self$generateScript()
      })
      gSignalConnect(inherit.lineheight.button, "toggled", function(obj){
        lineheight.inherit <<- obj$getActive()
        lineheight.button$setSensitive( !lineheight.inherit )
        .self$generateScript()
      })
      
      inherit.family.button$setActive(family.inherit)
      inherit.face.button$setActive(face.inherit)
      inherit.colour.button$setActive(colour.inherit)
      inherit.size.button$setActive(size.inherit)
      inherit.hjust.button$setActive(hjust.inherit)
      inherit.vjust.button$setActive(vjust.inherit)
      inherit.angle.button$setActive(angle.inherit)
      inherit.lineheight.button$setActive(lineheight.inherit)
      
      .self$generateScript()
    },
    
    buildCombo = function(contents, active=NULL){
      combo <- gtkComboBoxNewText()
      combo["width-request"] <- 1
      for(i in contents) combo$appendText(i)
      if (!is.null(active)) combo$setActive(active)
      return(combo)
    },
    
    generateScript = function(){
      if (name == "text" & family.inherit) {
        family.script <- sprintf("family = \"\"")
      } else {
        family.script <- sprintf("family = \"%s\"", if (family.inherit) NULL else family)        
      }
      face.script <- sprintf("face = \"%s\"", if (face.inherit) NULL else face)
      
      if (colour == "NA") {
        colour.script <- sprintf("colour = %s", if (colour.inherit) NULL else colour)
      } else {
        colour.script <- sprintf("colour = \"%s\"", if (colour.inherit) NULL else colour)        
      }
      
      if (size.rel) {
        size.script <- sprintf("size = rel(%s)", if (size.inherit) NULL else size)
      } else {
        size.script <- sprintf("size = %s", if (size.inherit) NULL else size)        
      }
      
      hjust.script      <- sprintf("hjust = %s"      , if (hjust.inherit     ) NULL else hjust     )
      vjust.script      <- sprintf("vjust = %s"      , if (vjust.inherit     ) NULL else vjust     )
      angle.script      <- sprintf("angle = %s"      , if (angle.inherit     ) NULL else angle     )
      lineheight.script <- sprintf("lineheight = %s" , if (lineheight.inherit) NULL else lineheight)
      script <<- paste(c(family.script, face.script, colour.script, size.script, hjust.script,
                         vjust.script, angle.script, lineheight.script), collapse=", ")
      if (blank) {
        script <<- sprintf("%s = element_blank()", name)
      } else if(nzchar(script)) {
        script <<- sprintf("%s = element_text(%s)", name, script)        
      } else {
        script <<- NULL
      }
      
    },
    
    setValue = function(value){
      if ("element_blank" %in% class(value)) {
        button.blank$setActive(FALSE)        
      } else {
        button.blank$setActive(TRUE)
        
        if (is.null(value$family)) {
          inherit.family.button$setActive(TRUE)          
        } else if (!nzchar(value$family)) {
          inherit.family.button$setActive(TRUE)
        } else {
          family.button$setFontName(value$family)
          inherit.family.button$setActive(FALSE)          
        }
        
        if (is.null(value$face)) {
          inherit.face.button$setActive(TRUE)          
        } else {
          face.combo$setActive( which(faces == value$face) - 1)
          inherit.face.button$setActive(FALSE)
        }
        
        if (is.null(value$colour)) {
          inherit.colour.button$setActive(TRUE)
        } else {
          color.widget$setColor(value$colour)
          inherit.colour.button$setActive(FALSE)
        }
        
        if (is.null(value$size)) {
          inherit.size.button$setActive(TRUE)          
        } else {
          size.button$setValue(value$size)
          if (class(value$size) == "rel") {
            size.rel.button$setActive(TRUE)
          } else {
            size.rel.button$setActive(FALSE)
          }
          inherit.size.button$setActive(FALSE)
        }
        
        if (is.null(value$hjust)) {
          inherit.hjust.button$setActive(TRUE)          
        } else {
          hjust.button$setValue(value$hjust)
          inherit.hjust.button$setActive(FALSE)
        }
        
        if (is.null(value$vjust)) {
          inherit.vjust.button$setActive(TRUE)          
        } else {
          vjust.button$setValue(value$vjust)
          inherit.vjust.button$setActive(FALSE)
        }
        
        if (is.null(value$angle)) {
          inherit.angle.button$setActive(TRUE)          
        } else {
          angle.button$setValue(value$angle)
          inherit.angle.button$setActive(FALSE)
        }
        
        if (is.null(value$lineheight)) {
          inherit.lineheight.button$setActive(TRUE)          
        } else {
          lineheight.button$setValue(value$lineheight)
          inherit.lineheight.button$setActive(FALSE)
        }
        
      }
    },
    
    reset = function(){
      family.inherit     <<- TRUE
      face.inherit       <<- TRUE
      colour.inherit     <<- TRUE
      size.inherit       <<- TRUE
      hjust.inherit      <<- TRUE
      vjust.inherit      <<- TRUE
      angle.inherit      <<- TRUE
      lineheight.inherit <<- TRUE
      inherit.family.button$setActive(family.inherit)
      inherit.face.button$setActive(face.inherit)
      inherit.colour.button$setActive(colour.inherit)
      inherit.size.button$setActive(size.inherit)
      inherit.hjust.button$setActive(hjust.inherit)
      inherit.vjust.button$setActive(vjust.inherit)
      inherit.angle.button$setActive(angle.inherit)
      inherit.lineheight.button$setActive(lineheight.inherit)
      button.blank$setActive(TRUE)
    }
  )
)
rzplot.element.text$accessors(c("main", "family", "face", "colour", "size",
                                "hjust", "vjust", "angle", "lineheight", "parent", "script"))
