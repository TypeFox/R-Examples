rzplot.element.rect <- 
setRefClass("RzPlotElementRect",
  fields = c("main", "name", "parent", "fill", "colour", "size", "size.rel", "linetype", "script", "blank",
             "fill.widget", "color.widget", "size.button", "size.rel.button", "linetype.combo",
             "inherit.from", "fill.inherit", "colour.inherit", "size.inherit", "linetype.inherit",
             "inherit.fill.button", "inherit.colour.button", "inherit.size.button", "inherit.linetype.button",
             "button.blank", "linetypes"),
  methods = list(
    initialize  = function(...) {
      initFields(...)
      fill.inherit     <<- FALSE
      colour.inherit   <<- FALSE
      size.inherit     <<- FALSE
      linetype.inherit <<- FALSE
      if (class(fill) == "uninitializedField") {
        fill <<- "#FFFFFF"
        fill.inherit <<- TRUE
      }
      if (class(colour) == "uninitializedField") {
        colour <<- "#FFFFFF"
        colour.inherit <<- TRUE
      }
      if (class(size) == "uninitializedField") {
        size <<- 0
        size.rel <<- FALSE
        size.inherit <<- TRUE
      }
      if (class(linetype) == "uninitializedField"){
        linetype <<- "blank"
        linetype.inherit <<- TRUE
      }
      if (class(blank) == "uninitializedField") {
        blank <<- FALSE
      }
      
      linetypes <<- c("blank", "solid", "dashed", "dotted",
                      "dotdash", "longdash", "twodash")
      
      fill.widget  <<- gtkColorSelectionWidgetNew(spacing=2, parent=parent)
      fill.widget$setColor(fill)
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
      linetype.combo <<- .self$buildCombo(linetypes, which(linetypes == linetype) - 1)
      
      table <- gtkTableNew()
      table$setBorderWidth(5)
      table$attach(gtkLabelNew(gettext("Fill"))     , 0, 1, 0, 1, "shrink", "shrink")
      table$attach(fill.widget                      , 1, 2, 0, 1, 5  , "shrink")
      table$attach(gtkLabelNew(gettext("Colour"))   , 0, 1, 1, 2, "shrink", "shrink")
      table$attach(color.widget                     , 1, 2, 1, 2, 5  , "shrink")
      table$attach(gtkLabelNew(gettext("Size"))     , 0, 1, 2, 3, "shrink", "shrink")
      table$attach(size.hbox                        , 1, 2, 2, 3, 5  , "shrink")
      table$attach(gtkLabelNew(gettext("Line\nType")), 0, 1, 3, 4, "shrink", "shrink")
      table$attach(linetype.combo                   , 1, 2, 3, 4, 5  , "shrink")
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
      
      gSignalConnect(fill.widget$getEntry(), "changed", function(obj){
        fill     <<- obj$getText()
        .self$generateScript()
      })
      gSignalConnect(color.widget$getEntry(), "changed", function(obj){
        colour   <<- obj$getText()
        .self$generateScript()
      })
      gSignalConnect(size.button, "changed", function(obj){
        size     <<- obj$getValue()
        .self$generateScript()
      })
      gSignalConnect(size.rel.button, "toggled", function(obj){
        size.rel <<- obj$getActive()
        .self$generateScript()
      })
      gSignalConnect(linetype.combo, "changed", function(obj){
        linetype <<- obj$getActiveText()
        .self$generateScript()
      })
      
      tooltip.text <- character()
      if (is.null(inherit.from)) {
        tooltip.text <- gettext("Inherit from the current theme")
      } else {
        tooltip.text <- gettextf("Inherit from <span font_style='italic' font_weight='bold' font_size='large'>%s</span>", inherit.from)
      }
      
      inherit.fill.button     <<- gtkToggleButtonNew()
      inherit.colour.button   <<- gtkToggleButtonNew()
      inherit.size.button     <<- gtkToggleButtonNew()
      inherit.linetype.button <<- gtkToggleButtonNew()
      inherit.fill.button$setImage(gtkImageNewFromStock(GTK_STOCK_GO_DOWN, GtkIconSize["button"]))
      inherit.colour.button$setImage(gtkImageNewFromStock(GTK_STOCK_GO_DOWN, GtkIconSize["button"]))
      inherit.size.button$setImage(gtkImageNewFromStock(GTK_STOCK_GO_DOWN, GtkIconSize["button"]))
      inherit.linetype.button$setImage(gtkImageNewFromStock(GTK_STOCK_GO_DOWN, GtkIconSize["button"]))
      inherit.fill.button["tooltip-markup"]     <<- tooltip.text
      inherit.colour.button["tooltip-markup"]   <<- tooltip.text
      inherit.size.button["tooltip-markup"]     <<- tooltip.text
      inherit.linetype.button["tooltip-markup"] <<- tooltip.text

      table$attach(inherit.fill.button    , 2, 3, 0, 1, "shrink", "shrink")
      table$attach(inherit.colour.button  , 2, 3, 1, 2, "shrink", "shrink")
      table$attach(inherit.size.button    , 2, 3, 2, 3, "shrink", "shrink")
      table$attach(inherit.linetype.button, 2, 3, 3, 4, "shrink", "shrink")
      
      gSignalConnect(inherit.fill.button, "toggled", function(obj){
        fill.inherit <<- obj$getActive()
        fill.widget$setSensitive( !fill.inherit )
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
      gSignalConnect(inherit.linetype.button, "toggled", function(obj){
        linetype.inherit <<- obj$getActive()
        linetype.combo$setSensitive( !linetype.inherit )
        .self$generateScript()
      })
      
      inherit.fill.button$setActive(fill.inherit)
      inherit.colour.button$setActive(colour.inherit)
      inherit.size.button$setActive(size.inherit)
      inherit.linetype.button$setActive(linetype.inherit)
      
      
      .self$generateScript()
    },
    
    buildCombo = function(contents, active=NULL){
      combo <- gtkComboBoxNewText()
      for(i in contents) combo$appendText(i)
      if (!is.null(active)) combo$setActive(active)
      return(combo)
    },
    
    generateScript = function(){
      if (fill == "NA") {
        fill.script <- sprintf("fill = %s"    , if (fill.inherit) NULL else fill)        
      } else {
        fill.script <- sprintf("fill = \"%s\"", if (fill.inherit) NULL else fill)        
      }
      
      if (colour == "NA") {
        colour.script <- sprintf("colour = %s"  , if (colour.inherit) NULL else colour)
      } else {
        colour.script <- sprintf("colour = \"%s\""  , if (colour.inherit) NULL else colour)        
      }
      
      if (size.rel) {
        size.script <- sprintf("size = rel(%s)", if (size.inherit) NULL else size)
      } else {
        size.script <- sprintf("size = %s", if (size.inherit) NULL else size)        
      }
      
      linetype.script <- sprintf("linetype = \"%s\"", if (linetype.inherit) NULL else linetype)
      script <<- paste(c(fill.script, colour.script, size.script, linetype.script), collapse=", ")
      if (blank) {
        script <<- sprintf("%s = element_blank()", name)
      } else if(nzchar(script)) {
        script <<- sprintf("%s = element_rect(%s)", name, script)        
      } else {
        script <<- NULL
      }
    },
    
    setValue = function(value){
      if ("element_blank" %in% class(value)) {
        button.blank$setActive(FALSE)        
      } else {
        button.blank$setActive(TRUE)
        
        if (is.null(value$fill)) {
          inherit.fill.button$setActive(TRUE)          
        } else {
          fill.widget$setColor(value$fill)
          inherit.fill.button$setActive(FALSE)          
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
        
        if (is.null(value$linetype)) {
          inherit.linetype.button$setActive(TRUE)          
        } else {
          if (is.numeric(value$linetype)) {
            linetype.combo$setActive(value$linetype)
          } else {
            linetype.combo$setActive( which(linetypes == value$linetype) - 1)
          }
          inherit.linetype.button$setActive(FALSE)
        }
      }
    },
    
    reset = function(){
      fill.inherit     <<- TRUE
      colour.inherit   <<- TRUE
      size.inherit     <<- TRUE
      linetype.inherit <<- TRUE
      inherit.fill.button$setActive(fill.inherit)
      inherit.colour.button$setActive(colour.inherit)
      inherit.size.button$setActive(size.inherit)
      inherit.linetype.button$setActive(linetype.inherit)
      button.blank$setActive(TRUE)
    }
  )
)
rzplot.element.rect$accessors(c("main", "fill", "colour", "size", "linetype", "script"))
