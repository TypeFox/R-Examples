VVDescriptives <-
setRefClass("RzVVDescriptives",
  fields = c("main"),
  contains = "RzVVCore",
  methods = list(
    initialize            = function(...) {
      initFields(...)      
      callSuper(...)
      
      main <<- gtkHPaned()
      rzAnalysisStat <- new("RzAnalysisStat", data=data, liststore=liststore, vvcore=.self)
      main$setPosition(400)
      main$pack1(widget, resize=TRUE)
      main$pack2(rzAnalysisStat$getMain(), resize=FALSE)
      
    },
    
    setAccel = function(accel.group){
#      button.execute$setAccelPath("<Rz-Menu>/View/Quick Editor View/Execute", accel.group)
    }
  )
)
VVDescriptives$accessors("main")
