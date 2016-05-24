VVPlot <-
setRefClass("RzVVPlot",
  fields = c("main"),
  contains ="RzVVCore",
  methods = list(
    initialize            = function(...) {
      initFields(...)
      callSuper(...)
      
      main <<- gtkHPaned()
      rzPlot <- new("RzPlot", vvcore=.self)
      rzPlot$setModel(liststore)
      rzPlot$setData(data)
      
      main$setPosition(400)
      main$pack1(widget, resize=TRUE)
      main$pack2(rzPlot$getMain(), resize=FALSE)
      rzPlot$construct()
      
    }
    
  )
)
VVPlot$accessors("main")
