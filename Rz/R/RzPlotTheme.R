rzplot.theme <- 
setRefClass("RzPlotTheme",
  fields = c("main",
             "line","rect","text","title","axis.title","axis.title.x","axis.title.y","axis.text",
             "axis.text.x","axis.text.y","axis.ticks","axis.ticks.x","axis.ticks.y",
             "axis.ticks.length","axis.ticks.margin","axis.line","axis.line.x","axis.line.y",
             "legend.background","legend.margin","legend.key","legend.key.size",
             "legend.key.height","legend.key.width","legend.text","legend.text.align",
             "legend.title","legend.title.align","legend.position","legend.direction",
             "legend.justification","legend.box","panel.background","panel.border","panel.margin",
             "panel.grid","panel.grid.major","panel.grid.minor","panel.grid.major.x",
             "panel.grid.major.y","panel.grid.minor.x","panel.grid.minor.y","plot.background",
             "plot.title","plot.margin","strip.background","strip.text","strip.text.x",
             "strip.text.y", "preview.path", "combo.load", "widget.list", "image.preview"),
  methods = list(
    initialize  = function(...) {
      initFields(...)
      if (! exists("theme_rz", envir=.GlobalEnv)) {
        rzTools$sync("theme_rz", theme_grey)        
      }
      
      main <<- gtkWindowNew(show=FALSE)
      main$setDefaultSize(1080, 640)
      main["title"] <<- gettext("Theme Editor")
      
      preview.path <<- file.path(tempdir(), "RzPlotThemePreview.png")
      
      gSignalConnect(main, "delete-event", function(obj, ...) obj$hideOnDelete())
      
      #       line   all line elements (element_line)
      #       rect   all rectangluar elements (element_rect)
      #       text	 all text elements (element_text)
      #       title	 all title elements: plot, axes, legends (element_text; inherits from text)
      line  <<- new("RzPlotElementLine", name="line",  parent=main, inherit.from=NULL)
      rect  <<- new("RzPlotElementRect", name="rect",  parent=main, inherit.from=NULL)
      text  <<- new("RzPlotElementText", name="text",  parent=main, inherit.from=NULL)
      title <<- new("RzPlotElementText", name="title", parent=main, inherit.from=NULL)
      
      #       axis.title  label of axes (element_text; inherits from text)
      #       axis.text   tick labels along axes (element_text; inherits from text)
      #       axis.ticks   tick marks along axes (element_line; inherits from line)      
      #       axis.ticks.length   length of tick marks (unit)
      #       axis.ticks.margin	 space between tick mark and tick label (unit)
      #       axis.line   lines along axes (element_line; inherits from line)
      axis.title        <<- new("RzPlotElementText", name="axis.title", parent=main, inherit.from="text")
      axis.text         <<- new("RzPlotElementText", name="axis.text", parent=main, inherit.from="text")
      axis.ticks        <<- new("RzPlotElementLine", name="axis.ticks", parent=main, inherit.from="line")
      axis.ticks.length <<- gtkPlotUnitNew(name="axis.ticks.length")
      axis.ticks.margin <<- gtkPlotUnitNew(name="axis.ticks.margin")
      axis.line         <<- new("RzPlotElementLine", name="axis.line", parent=main, inherit.from="line")
      
      #       axis.title.x	 x axis label (element_text; inherits from axis.title)
      #       axis.text.x   x axis tick labels (element_text; inherits from axis.text)
      #       axis.ticks.x   x axis tick marks (element_line; inherits from axis.ticks)
      #       axis.line.x   line along x axis (element_line; inherits from axis.line)
      axis.title.x <<- new("RzPlotElementText", name="axis.title.x", parent=main, inherit.from="axis.title")
      axis.text.x  <<- new("RzPlotElementText", name="axis.text.x", parent=main, inherit.from="axis.text")
      axis.ticks.x <<- new("RzPlotElementLine", name="axis.ticks.x", parent=main, inherit.from="axis.ticks")
      axis.line.x  <<- new("RzPlotElementLine", name="axis.line.x", parent=main, inherit.from="axis.line")
      
      #       axis.title.y	 y axis label (element_text; inherits from axis.title)
      #       axis.text.y	 y axis tick labels (element_text; inherits from axis.text)
      #       axis.ticks.y	 y axis tick marks (element_line; inherits from axis.ticks)
      #       axis.line.y	 line along y axis (element_line; inherits from axis.line)
      axis.title.y <<- new("RzPlotElementText", name="axis.title.y", parent=main, inherit.from="axis.title")
      axis.text.y  <<- new("RzPlotElementText", name="axis.text.y", parent=main, inherit.from="axis.text")
      axis.ticks.y <<- new("RzPlotElementLine", name="axis.ticks.y", parent=main, inherit.from="axis.ticks")
      axis.line.y  <<- new("RzPlotElementLine", name="axis.line.y", parent=main, inherit.from="axis.line")     
      
      #       legend.background  background of legend (element_rect; inherits from rect)
      #       legend.margin	 extra space added around legend (unit)
      #       legend.key	 background underneath legend keys (element_rect; inherits from rect)
      #       legend.key.size	 size of legend keys (unit; inherits from legend.key.size)
      #       legend.key.height	 key background height (unit; inherits from legend.key.size)
      #       legend.key.width	 key background width (unit; inherits from legend.key.size)
      #       legend.text	 legend item labels (element_text; inherits from text)
      #       legend.text.align	 alignment of legend labels (number from 0 (left) to 1 (right))
      #       legend.title	 title of legend (element_text; inherits from title)
      #       legend.title.align	 alignment of legend title (number from 0 (left) to 1 (right))
      #       legend.position	 the position of legends. ("left", "right", "bottom", "top", or two-element numeric vector)
      #       legend.direction	 layout of items in legends ("horizontal" or "vertical")
      #       legend.justification	 anchor point for positioning legend inside plot ("center" or two-element numeric vector)
      #       legend.box	 arrangement of multiple legends ("horizontal" or "vertical")
      legend.background    <<- new("RzPlotElementRect", name="legend.background",  parent=main, inherit.from="rect")
      legend.margin        <<- gtkPlotUnitNew(name="legend.margin")
      legend.key           <<- new("RzPlotElementRect", name="legend.key",  parent=main, inherit.from="rect")
      legend.key.size      <<- gtkPlotUnitNew(name="legend.key.size")
      legend.key.height    <<- gtkPlotUnitNew(name="legend.key.height", inherit.from="legend.key.size")
      legend.key.width     <<- gtkPlotUnitNew(name="legend.key.width", inherit.from="legend.key.size")
      legend.text          <<- new("RzPlotElementText", name="legend.text", parent=main, inherit.from="text")
      legend.text.align    <<- gtkPlotAlignNew(name="legend.text.align")
      legend.title         <<- new("RzPlotElementText", name="legend.title", parent=main, inherit.from="title")
      legend.title.align   <<- gtkPlotAlignNew(name="legend.title.align")
      legend.position      <<- gtkPlotPositionNew(name="legend.position")
      legend.direction     <<- gtkPlotDirectionNew(name="legend.direction")
      legend.justification <<- gtkPlotPositionNew(name="legend.justification")
      legend.box           <<- gtkPlotDirectionNew(name="legend.box")
      
      #       panel.background  background of plotting area, drawn underneath plot (element_rect; inherits from rect)
      #       panel.border	 border around plotting area, drawn on top of plot so that it covers tick marks and grid lines. This should be used with fill=NA (element_rect; inherits from rect)
      #       panel.margin	 margin around facet panels (unit)
      #       panel.grid	 grid lines (element_line; inherits from line)
      #       panel.grid.major	 major grid lines (element_line; inherits from panel.grid)
      #       panel.grid.minor	 minor grid lines (element_line; inherits from panel.grid)
      #       panel.grid.major.x	 vertical major grid lines (element_line; inherits from panel.grid.major)
      #       panel.grid.major.y	 horizontal major grid lines (element_line; inherits from panel.grid.major)
      #       panel.grid.minor.x	 vertical minor grid lines (element_line; inherits from panel.grid.minor)
      #       panel.grid.minor.y	 horizontal minor grid lines (element_line; inherits from panel.grid.minor)
      panel.background   <<- new("RzPlotElementRect", name="panel.background",  parent=main, inherit.from="rect")
      panel.border       <<- new("RzPlotElementRect", name="panel.border",  parent=main, inherit.from="rect")
      panel.margin       <<- gtkPlotUnitNew(name="panel.margin")
      panel.grid         <<- new("RzPlotElementLine", name="panel.grid", parent=main, inherit.from="line")
      panel.grid.major   <<- new("RzPlotElementLine", name="panel.grid.major", parent=main, inherit.from="panel.grid")
      panel.grid.minor   <<- new("RzPlotElementLine", name="panel.grid.minor", parent=main, inherit.from="panel.grid")
      panel.grid.major.x <<- new("RzPlotElementLine", name="panel.grid.major.x", parent=main, inherit.from="panel.grid.major")
      panel.grid.major.y <<- new("RzPlotElementLine", name="panel.grid.major.y", parent=main, inherit.from="panel.grid.major")
      panel.grid.minor.x <<- new("RzPlotElementLine", name="panel.grid.minor.x", parent=main, inherit.from="panel.grid.minor")
      panel.grid.minor.y <<- new("RzPlotElementLine", name="panel.grid.minor.y", parent=main, inherit.from="panel.grid.minor")
      
      #       plot.background  background of the entire plot (element_rect; inherits from rect)
      #       plot.title	 plot title (text appearance) (element_text; inherits from title)
      #       plot.margin	 margin around entire plot (unit with the sizes of the top, right, bottom, and left margins)
      plot.background <<- new("RzPlotElementRect", name="plot.background",  parent=main, inherit.from="rect")
      plot.title      <<- new("RzPlotElementText", name="plot.title", parent=main, inherit.from="title")
      plot.margin     <<- gtkPlotUnit2New(name="plot.margin")
      
      #       strip.background  background of facet labels (element_rect; inherits from rect)
      #       strip.text	 facet labels (element_text; inherits from text)
      #       strip.text.x	 facet labels along horizontal direction (element_text; inherits from strip.text)
      #       strip.text.y	 facet labels along vertical direction (element_text; inherits from strip.text)
      strip.background <<- new("RzPlotElementRect", name="strip.background",  parent=main, inherit.from="rect")
      strip.text       <<- new("RzPlotElementText", name="strip.text", parent=main, inherit.from="text")
      strip.text.x     <<- new("RzPlotElementText", name="strip.text.x", parent=main, inherit.from="strip.text")
      strip.text.y     <<- new("RzPlotElementText", name="strip.text.y", parent=main, inherit.from="strip.text")
      
      widget.list <<- list(line, rect, text, title, axis.title, axis.title.x, axis.title.y, axis.text, 
                           axis.text.x, axis.text.y, axis.ticks, axis.ticks.x, axis.ticks.y, 
                           axis.ticks.length, axis.ticks.margin, axis.line, axis.line.x, axis.line.y, 
                           legend.background, legend.margin, legend.key, legend.key.size, legend.key.height, 
                           legend.key.width, legend.text, legend.text.align, legend.title, 
                           legend.title.align, legend.position, legend.direction, legend.justification, 
                           legend.box, panel.background, panel.border, panel.margin, panel.grid, 
                           panel.grid.major, panel.grid.minor, panel.grid.major.x, panel.grid.major.y, 
                           panel.grid.minor.x, panel.grid.minor.y, plot.background, plot.title, plot.margin, 
                           strip.background, strip.text, strip.text.x, strip.text.y)
      names(widget.list) <<- c("line","rect","text","title","axis.title","axis.title.x","axis.title.y","axis.text",
                               "axis.text.x","axis.text.y","axis.ticks","axis.ticks.x","axis.ticks.y",
                               "axis.ticks.length","axis.ticks.margin","axis.line","axis.line.x","axis.line.y",
                               "legend.background","legend.margin","legend.key","legend.key.size",
                               "legend.key.height","legend.key.width","legend.text","legend.text.align",
                               "legend.title","legend.title.align","legend.position","legend.direction",
                               "legend.justification","legend.box","panel.background","panel.border","panel.margin",
                               "panel.grid","panel.grid.major","panel.grid.minor","panel.grid.major.x",
                               "panel.grid.major.y","panel.grid.minor.x","panel.grid.minor.y","plot.background",
                               "plot.title","plot.margin","strip.background","strip.text","strip.text.x",
                               "strip.text.y")

      table.general <- gtkTableNew()
      table.general$setColSpacings(2)
      table.general$setRowSpacings(5)
      table.general$attach(line$getMain(), 0, 1, 0, 1)
      table.general$attach(rect$getMain(), 1, 2, 0, 1)
      table.general$attach(text$getMain(), 0, 1, 1, 2)
      table.general$attach(title$getMain(),1, 2, 1, 2)
      sw.general <- gtkScrolledWindowWithViewportNew()
      sw.general$add(table.general)
      
      table.axis <- gtkTableNew()
      table.axis$setColSpacings(2)
      table.axis$setRowSpacings(5)
      table.axis$attach(axis.title$getMain(),0, 1, 0, 1)
      table.axis$attach(axis.text$getMain() ,1, 2, 0, 1)
      table.axis$attach(axis.line$getMain(), 0, 1, 1, 2)
      table.axis$attach(axis.ticks$getMain(),1, 2, 1, 2)
      table.axis$attach(axis.ticks.length,   0, 1, 2, 3)
      table.axis$attach(axis.ticks.margin,   1, 2, 2, 3)
      sw.axis <- gtkScrolledWindowWithViewportNew()
      sw.axis$add(table.axis)
      
      table.axis.x <- gtkTableNew()
      table.axis.x$setColSpacings(2)
      table.axis.x$setRowSpacings(5)
      table.axis.x$attach(axis.title.x$getMain(), 0, 1, 0, 1)
      table.axis.x$attach(axis.text.x$getMain(),  1, 2, 0, 1)
      table.axis.x$attach(axis.line.x$getMain(),  0, 1, 1, 2)
      table.axis.x$attach(axis.ticks.x$getMain(), 1, 2, 1, 2)
      sw.axis.x <- gtkScrolledWindowWithViewportNew()
      sw.axis.x$add(table.axis.x)
      
      table.axis.y <- gtkTableNew()
      table.axis.y$setColSpacings(2)
      table.axis.y$setRowSpacings(5)
      table.axis.y$attach(axis.title.y$getMain(), 0, 1, 0, 1)
      table.axis.y$attach(axis.text.y$getMain(),  1, 2, 0, 1)
      table.axis.y$attach(axis.line.y$getMain(),  0, 1, 1, 2)
      table.axis.y$attach(axis.ticks.y$getMain(), 1, 2, 1, 2)
      sw.axis.y <- gtkScrolledWindowWithViewportNew()
      sw.axis.y$add(table.axis.y)
            
      table.legend <- gtkTableNew()
      table.legend$setColSpacings(2)
      table.legend$setRowSpacings(5)
      table.legend$attach(legend.text$getMain(),       0, 1, 0, 1)
      table.legend$attach(legend.title$getMain(),      1, 2, 0, 1)
      table.legend$attach(legend.margin,               0, 1, 1, 2)
      table.legend$attach(legend.text.align,           0, 1, 2, 3)
      table.legend$attach(legend.title.align,          0, 1, 3, 4)
      table.legend$attach(legend.background$getMain(), 1, 2, 1, 4)
      table.legend$attach(legend.key$getMain(),        0, 1, 4, 7)
      table.legend$attach(legend.key.size,             1, 2, 4, 5)
      table.legend$attach(legend.key.height,           1, 2, 5, 6)
      table.legend$attach(legend.key.width,            1, 2, 6, 7)
      table.legend$attach(legend.position,             0, 1, 7, 8)
      table.legend$attach(legend.justification,        1, 2, 7, 8)
      table.legend$attach(legend.direction,            0, 1, 8, 9)
      table.legend$attach(legend.box,                  1, 2, 8, 9)      
      sw.legend <- gtkScrolledWindowWithViewportNew()
      sw.legend$add(table.legend)      
      
      table.panel <- gtkTableNew()
      table.panel$setColSpacings(2)
      table.panel$setRowSpacings(5)
      table.panel$attach(panel.background$getMain(),   0, 1, 0, 1)
      table.panel$attach(panel.border$getMain(),       1, 2, 0, 1)
      table.panel$attach(panel.margin,                 0, 1, 1, 2)
      table.panel$attach(panel.grid$getMain(),         1, 2, 1, 2)
      table.panel$attach(panel.grid.major$getMain(),   0, 1, 2, 3)
      table.panel$attach(panel.grid.minor$getMain(),   1, 2, 2, 3)
      table.panel$attach(panel.grid.major.x$getMain(), 0, 1, 3, 4)
      table.panel$attach(panel.grid.major.y$getMain(), 1, 2, 3, 4)
      table.panel$attach(panel.grid.minor.x$getMain(), 0, 1, 4, 5)
      table.panel$attach(panel.grid.minor.y$getMain(), 1, 2, 4, 5)
      sw.panel <- gtkScrolledWindowWithViewportNew()
      sw.panel$add(table.panel)
      
      table.plot <- gtkTableNew()
      table.plot$setColSpacings(2)
      table.plot$setRowSpacings(5)
      table.plot$attach(plot.title$getMain(),      0, 1, 0, 2)
      table.plot$attach(plot.background$getMain(), 1, 2, 0, 1)
      table.plot$attach(plot.margin,               1, 2, 1, 2)
      sw.plot <- gtkScrolledWindowWithViewportNew()
      sw.plot$add(table.plot)
      
      table.strip <- gtkTableNew()
      table.strip$setColSpacings(2)
      table.strip$setRowSpacings(5)
      table.strip$attach(strip.background$getMain(), 0, 1, 0, 1)
      table.strip$attach(strip.text$getMain(),       1, 2, 0, 1)
      table.strip$attach(strip.text.x$getMain(),     0, 1, 1, 2)
      table.strip$attach(strip.text.y$getMain(),     1, 2, 1, 2)
      sw.strip <- gtkScrolledWindowWithViewportNew()
      sw.strip$add(table.strip)
      
      note <- gtkNotebookNew()
      note$setScrollable(TRUE)
      note$appendPage(sw.general, tab.label=gtkLabelNew(gettext("General")))
      note$appendPage(sw.axis,    tab.label=gtkLabelNew(gettext("Axis")))
      note$appendPage(sw.axis.x,  tab.label=gtkLabelNew(gettext("X Axis")))
      note$appendPage(sw.axis.y,  tab.label=gtkLabelNew(gettext("Y Axis")))
      note$appendPage(sw.legend,  tab.label=gtkLabelNew(gettext("Legend")))
      note$appendPage(sw.panel,   tab.label=gtkLabelNew(gettext("Panel")))
      note$appendPage(sw.plot,    tab.label=gtkLabelNew(gettext("Plot")))
      note$appendPage(sw.strip,   tab.label=gtkLabelNew(gettext("Strip")))
      
      
      # Controll
      
      combo.load <<- gtkComboBoxNewText()
      combo.load["width-request"] <<- 120
      combo.load$appendText("theme_rz")
      combo.load$appendText("theme_grey")
      combo.load$appendText("theme_bw")
      if (packageVersion("ggplot2") >= "0.9.3") {
        combo.load$appendText("theme_minimal")
        combo.load$appendText("theme_classic")        
      }
      combo.load$setActive(0)

      button.load     <- gtkButtonNewWithLabel("Load")
      button.open     <- gtkButtonNewWithLabel("Open...")
      button.save     <- gtkButtonNewWithLabel("Save...")
      button.reset    <- gtkButtonNewWithLabel("Reset")
      button.refresh1 <- gtkButtonNewWithLabel("Preview 1")
      button.refresh2 <- gtkButtonNewWithLabel("Preview 2")
      button.refresh3 <- gtkButtonNewWithLabel("Preview 3")
      button.refresh4 <- gtkButtonNewWithLabel("Preview 4")
      button.output   <- gtkButtonNewWithLabel("Output")
      button.apply    <- gtkButtonNewWithLabel("Apply")
      button.ok       <- gtkButtonNewWithLabel("Ok")
      
      sizegroup <- gtkSizeGroupNew()
      sizegroup$addWidget(button.load)
      sizegroup$addWidget(button.open)
      sizegroup$addWidget(button.save)
      sizegroup$addWidget(button.reset)
      sizegroup$addWidget(button.output)
      sizegroup$addWidget(button.apply)
      sizegroup$addWidget(button.ok)
            
      hbox.preview <- gtkHBoxNew(spacing=2)
      hbox.preview$packEnd(button.refresh4, expand=FALSE)
      hbox.preview$packEnd(button.refresh3, expand=FALSE)
      hbox.preview$packEnd(button.refresh2, expand=FALSE)
      hbox.preview$packEnd(button.refresh1, expand=FALSE)
      
      df  <- data.frame(x = 1:1000, y=1:1000)
      plot.preview <- ggplot(data=df, aes(x, y)) + geom_blank() +
        scale_x_continuous(breaks=c(0, 250, 500, 750, 1000),
                           labels=c(0, 250, "axis.text.x", 750, 1000)) +
        scale_y_continuous(breaks=c(0, 250, 500, 750, 1000),
                           labels=c(0, 250, "axis.text.y", 750, 1000)) +
        labs(title="plot.title", x="axis.title.x", y="axis.title.y") +
        annotate("text", x = 500, y = 500, label = "panel")
      suppressMessages(ggsave(preview.path, plot.preview, width=7, height=7))
      pixbuf <- gdkPixbufNewFromFileAtSize(preview.path, width=480, height=-1)$retval
      image.preview <<- gtkImageNewFromPixbuf(pixbuf)
      sw.preview <- gtkScrolledWindowWithViewportNew()
      sw.preview$add(image.preview)
      
      vbox.preview <- gtkVBoxNew(spacing=5)
      vbox.preview$setBorderWidth(5)
      vbox.preview$packStart(sw.preview,   expand=TRUE)
      vbox.preview$packStart(hbox.preview, expand=FALSE)
      
      frame.preview <- gtkFrameNew(gettext("Theme Preview"))
      frame.preview$add(vbox.preview)
      
      hbox.save <- gtkHBoxNew(spacing=2)
      hbox.save$packStart(button.load, expand=FALSE)
      hbox.save$packStart(button.open, expand=FALSE)
      hbox.save$packStart(button.save, expand=FALSE)
      hbox.save$packEnd(button.apply,  expand=FALSE)
      hbox.save$packEnd(button.ok,     expand=FALSE)
      hbox.save$packEnd(button.output, expand=FALSE)
      
      vbox.save <- gtkVBoxNew()
      vbox.save$packStart(hbox.save, expand=TRUE, fill=FALSE)
      
      hbox.save2 <- gtkHBoxNew(spacing=2)
      hbox.save2$packStart(combo.load,  expand=FALSE)
      hbox.save2$packStart(vbox.save,   expand=TRUE)
      
      vbox.controll <- gtkVBoxNew(spacing=5)
      vbox.controll$packStart(frame.preview, expand=TRUE)
      vbox.controll$packEnd(hbox.save2,    expand=FALSE)
      
      vbox.note <- gtkVBoxNew()
      
      if(is.null(gtkCheckVersion(2, 20, 0))) {
        vbox.note$packStart(note)
        note$setActionWidget(button.reset, "end")        
      } else {
        hbox.button.reset <- gtkHBoxNew()
        hbox.button.reset$packEnd(button.reset, expand=FALSE)
        vbox.note$packStart(hbox.button.reset, expand=FALSE)
        vbox.note$packStart(note)
      }
      
      paned <- gtkHPanedNew()
      paned$pack1(vbox.note, resize=FALSE, shrink=FALSE)
      paned$pack2(vbox.controll, resize=TRUE, shrink=FALSE)
      paned$setPosition(550)
      paned$setBorderWidth(5)
      main$add(paned)
      
      gSignalConnect(button.refresh1, "clicked", .self$onRefresh, 1)
      gSignalConnect(button.refresh2, "clicked", .self$onRefresh, 2)
      gSignalConnect(button.refresh3, "clicked", .self$onRefresh, 3)
      gSignalConnect(button.refresh4, "clicked", .self$onRefresh, 4)
      gSignalConnect(button.load,     "clicked", .self$onLoad)
      gSignalConnect(button.open,     "clicked", .self$onOpen)
      gSignalConnect(button.save,     "clicked", .self$onSave)
      gSignalConnect(button.reset ,   "clicked", .self$onReset)
      gSignalConnect(button.output,   "clicked", .self$onOutput)
      gSignalConnect(button.apply,    "clicked", .self$onApply, FALSE)
      gSignalConnect(button.ok,       "clicked", .self$onApply, TRUE)
      
      .self$onLoad()
    },
    
    show = function(){
      main$show()
      main$present()
    },
    
    onRefresh = function(object, data){
      main$setKeepAbove(TRUE)
      script <- .self$getScript()
      script[[2]] <- paste(script[[2]], collapse=", ")

      if (!is.null(script[[1]])) {
        eval(parse(text = script[[1]]))
      }

      theme.tmp <- eval(parse(text = sprintf("theme(%s, complete=TRUE)", script[[2]])))
      if (data==1) {
        df  <- data.frame(x = 1:1000, y=1:1000)
        plot.preview <- ggplot(data=df, aes(x, y)) + geom_blank() +
          scale_x_continuous(breaks=c(0, 250, 500, 750, 1000),
                             labels=c(0, 250, "axis.text.x", 750, 1000)) +
          scale_y_continuous(breaks=c(0, 250, 500, 750, 1000),
                             labels=c(0, 250, "axis.text.y", 750, 1000)) +
          labs(title="plot.title", x="axis.title.x", y="axis.title.y") +
          annotate("text", x = 500, y = 500, label = "panel") + theme.tmp
        suppressMessages(ggsave(preview.path, plot.preview, width=7, height=7))
        pixbuf <- gdkPixbufNewFromFileAtSize(preview.path, width=480, height=-1)$retval
      } else if (data==2) {
        df <- data.frame(x  = 1:100,
                         y  = 1:100,
                         sx = rep(c("strip.text.x1", "strip.text.x2"), 50),
                         sy = rep(c("strip.text.y1", "strip.text.y2"), 50),
                         l  = rep(c("legend.text1", "legend.text2"), 50))
        plot.preview <- ggplot(data=df, aes(x, y, color=l, shape=l)) + geom_point(size=0) +
          facet_grid(sy~sx) + 
          guides(color=guide_legend(title="legend.title1", order=1),
                 shape=guide_legend(title="legend.title2", order=0)) +
          scale_x_continuous(breaks=c(15, 50, 85),
                             labels=c("axis.text.x1", "axis.text.x2", "axis.text.x3")) +
          scale_y_continuous(breaks=c(15, 50, 85),
                             labels=c("axis.text.y1", "axis.text.y2", "axis.text.y3")) +
          labs(title="plot.title", x="axis.title.x", y="axis.title.y") +
          annotate("text", x = 50, y = 50, label = "panel") + theme.tmp
        suppressMessages(ggsave(preview.path, plot.preview, width=8, height=6))
        pixbuf <- gdkPixbufNewFromFileAtSize(preview.path, width=600, height=-1)$retval
      } else if (data==3) {
        data(iris, envir=environment())
        plot.preview <- ggplot(iris, aes(Sepal.Length, Petal.Length)) + 
          geom_point() + 
          labs(title="plot.title", x="axis.title.x", y="axis.title.y") + theme.tmp
        suppressMessages(ggsave(preview.path, plot.preview, width=7, height=7))
        pixbuf <- gdkPixbufNewFromFileAtSize(preview.path, width=480, height=-1)$retval
      } else {
        df <- subset(diamonds, subset=
                       cut %in% c("Good", "Very Good") &
                       color %in% c("D","F") &
                       clarity %in% c("VVS2", "IF"),
                     select = c(carat:clarity, price))
        df$cut     <- df$cut[drop=TRUE]
        df$color   <- df$color[drop=TRUE]
        df$clarity <- df$clarity[drop=TRUE]
        df$price   <- df$price / 10000
        levels(df$cut)     <- c("strip.text.x1", "strip.text.x2")
        levels(df$color)   <- c("legend.text1", "legend.text2")
        levels(df$clarity) <- c("strip.text.y1", "strip.text.y2")
        colnames(df) <- c("axis.title.y", "cut", "legend.title", "clarity", "axis.title.x")
        plot.preview <- ggplot(df, aes(axis.title.x, axis.title.y, color=legend.title, shape=legend.title)) +
          geom_point() + facet_grid(clarity ~ cut) +
          guides(color=guide_legend(title="legend.title1", order=1),
                 shape=guide_legend(title="legend.title2", order=0)) +
          labs(title="plot.title", x="axis.title.x", y="axis.title.y") + theme.tmp
        suppressMessages(ggsave(preview.path, plot.preview, width=8, height=6))
        pixbuf <- gdkPixbufNewFromFileAtSize(preview.path, width=600, height=-1)$retval
      }

      image.preview$setFromPixbuf(pixbuf)
      main$present()
      main$setKeepAbove(FALSE)
    },
    
    onLoad = function(...){
      .self$onReset()
      theme.name <- localize(combo.load$getActiveText())
      if (length(theme.name) == 0) return()
      
      theme.elements <- get(theme.name)()
      elements.names <- names(theme.elements)
      
      for (i in elements.names) {
        widget.list[[i]]$setValue(theme.elements[[i]])
      }
    },
    
    onOpen = function(...){
      dialog <- gtkFileChooserDialogNew(title=gettext("Open"), parent=main,
                                        action=GtkFileChooserAction["open"],
                                        "gtk-save", GtkResponseType["accept"],
                                        "gtk-cancel", GtkResponseType["cancel"], 
                                        show=FALSE)
      filter <- gtkFileFilterNew()
      filter$setName("R Script (*.R)")
      filter$addPattern("*.R")
      dialog$addFilter(filter)
      
      gSignalConnect(dialog, "response", function(object, response.id, ...){
        if ( response.id == GtkResponseType["accept"] ) {
          filename <- localize(dialog$getFilename())
          env.tmp <- new.env()
          source(filename, local=env.tmp)
          theme.fun.name <- ls(envir=env.tmp)
          theme.tmp <- try(env.tmp[[theme.fun.name]](), silent=TRUE)
          if (all(is.theme(theme.tmp))) {
            source(filename)
            combo.load$appendText(theme.fun.name)
            model <- combo.load$getModel()
            iter  <- model$getIterFirst()
            active <- NULL
            while (iter$retval) {
              active      <- model$getStringFromIter(iter$iter)
              iter$retval <- model$iterNext(iter$iter)
            }
            combo.load$setActive(active)
            .self$onLoad()
          }
        } else {
          return()
        }
      })
      
      dialog$run()
      dialog$hide()
    },
    
    onReset = function(...){
      lapply(widget.list, function(x) x$reset())
    },
    
    getScript = function(...){
      fonts.script <- NULL
      script <- lapply(widget.list, function(x) x$getScript())
      script <- unlist(script)
      fonts <- grep("family = ", script, value=TRUE)
      fonts <- sub("(^.*family = \\\")([^,]*)(\\\".*$)", "\\2", fonts)
      fonts <- unique(fonts)
      fonts <- fonts[nzchar(fonts)]
      if (length(fonts) != 0) {
        fonts.script <- fontsRegisterScript(fonts)        
      }
      return(list(fonts.script, script))
    },
    
    onOutput = function(...){
      script <- .self$getScript()
      if (!is.null(script[[1]])) {
        cat(script[[1]], fill=TRUE)        
      }
      cat("theme_rz <- function(){\n  theme(\n")
      cat(paste("   ", script[[2]]), "    complete = TRUE)", sep=",\n")
      cat("}", fill=TRUE)
    },
    
    onApply = function(object, data){
      script <- .self$getScript()
      script[[2]] <- paste(script[[2]], collapse=", ")
      eval(parse(text = sprintf("theme_rz <- function(){theme(%s, complete=TRUE)}", script[[2]])),
           envir=.GlobalEnv)
      
      if (!is.null(script[[1]])) {
        eval(parse(text = script[[1]]), envir=.GlobalEnv)
      }
      
      if (data) {
        main$hide()
      }
      
    },
    
    onSave = function(object){
      dialog <- gtkFileChooserDialogNew(title=gettext("Save"), parent=main,
                                        action=GtkFileChooserAction["save"],
                                        "gtk-save", GtkResponseType["accept"],
                                        "gtk-cancel", GtkResponseType["cancel"], 
                                        show=FALSE)
      filter <- gtkFileFilterNew()
      filter$setName("R Script (*.R)")
      filter$addPattern("*.R")
      dialog$addFilter(filter)
      
      label1 <- gtkLabelNew(gettext("Name of Theme :"))
      label2 <- gtkLabelNew(" theme_")
      entry  <- gtkEntryNew()
      entry["tooltip-text"] <- gettext("use filename if blank")
      hbox <- gtkHBoxNew()
      hbox$packStart(label1,expand=FALSE)
      hbox$packStart(label2,expand=FALSE)
      hbox$packStart(entry,expand=FALSE)
      dialog$setExtraWidget(hbox)
      
      gSignalConnect(dialog, "response", function(object, response.id, ...){
        if ( response.id == GtkResponseType["accept"] ) {
          filename <- localize(dialog$getFilename())
          if (! grepl("*\\.R$|*\\.r$", filename)) {
            filename <- sprintf("%s.R", filename)
          }
          
          if (file.exists(filename)){
            dialog2 <- gtkMessageDialogNew(dialog, "destroy-with-parent",
                                           GtkMessageType["question"],
                                           GtkButtonsType["ok-cancel"],
                                           gettext("This file already exists. Overwrite it?"))
            response <- dialog2$run()
            dialog2$hide()
            if (response!=GtkResponseType["ok"]){
              dialog$run()
              return()
            }
          }
          
          basename <- basename(filename)
          basename <- sub("\\.R", "", basename)
          themename <- entry$getText()
          if (nzchar(themename)) {
            basename <- themename
          }

          script <- .self$getScript()

          if (!is.null(script[[1]])) {
            cat(script[[1]], fill=TRUE, file=filename)
            cat(sprintf("theme_%s <- function(){\n  theme(\n", basename), file=filename, append=TRUE)
            cat(paste("   ", script[[2]]), "    complete = TRUE)", sep=",\n", file=filename, append=TRUE)
            cat("}\n", fill=TRUE, file=filename, append=TRUE)
          } else {
            cat(sprintf("theme_%s <- function(){\n  theme(\n", basename), file=filename)
            cat(paste("   ", script[[2]]), "    complete = TRUE)", sep=",\n", file=filename, append=TRUE)
            cat("}\n", fill=TRUE, file=filename, append=TRUE)            
          }
          
        } else {
          return()
        }
      })

      dialog$run()
      dialog$hide()
    }
  )
)
rzplot.theme$accessors("main")

#obj <- rzplot.theme$new()
#obj$getMain()$show()
