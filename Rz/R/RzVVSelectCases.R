selectCases <-
  setRefClass("RzVVSelectCases",
  fields = c("main", "vvcore", "data", "textView", "textBuffer", "button.update", "button.clear"),
  methods = list(
    initialize = function(...) {
      initFields(...)
      
      toggleButton  <-  gtkToggleButtonNewWithLabel(gettext("Enable Select Cases"))
      button.update <<- gtkButtonNewFromStock(GTK_STOCK_REFRESH)
      button.clear  <<- gtkButtonNewFromStock(GTK_STOCK_CLEAR)
      textBuffer    <<- gtkTextBufferNew()
      textBuffer$setText(data$getSubset.condition())
      textView      <<- gtkTextViewNewWithBuffer(textBuffer)
      textView$setLeftMargin(5)
      textView$setRightMargin(5)
      
      subset.on <- data$getSubset.on()
      toggleButton$setActive(subset.on)
      button.update$setSensitive(subset.on)
      button.clear$setSensitive(subset.on)
      textView$setSensitive(subset.on)
      
      scrolledWindow.textView <- gtkScrolledWindowNew()
      scrolledWindow.textView$setShadowType(GtkShadowType["in"])
      scrolledWindow.textView$setPolicy(GtkPolicyType["automatic"], GtkPolicyType["automatic"])
      scrolledWindow.textView$add(textView)
      scrolledWindow.textView$setSizeRequest(-1, 100)
      
      hbox.button1 <- gtkHBoxNew()
      hbox.button1$packStart(toggleButton, expand=FALSE)
      
      hbox.button2 <- gtkHBoxNew()
      hbox.button2$packEnd(button.update, expand=FALSE)
      hbox.button2$packStart(button.clear, expand=FALSE)

      vbox1 <- gtkVBoxNew(spacing=2)
      vbox1$packStart(hbox.button1, expand=FALSE)
      vbox1$packStart(scrolledWindow.textView, expand=FALSE)
      vbox1$packStart(hbox.button2, expand=FALSE)
      
      vbox <- gtkVBoxNew(spacing=4)
      vbox$setBorderWidth(2)
      vbox$packStart(vbox1, expand=FALSE)
      
      main <<- gtkScrolledWindowNew()
      main$setShadowType(GtkShadowType["none"])
      main$setPolicy(GtkPolicyType["automatic"], GtkPolicyType["automatic"])
      main$addWithViewport(vbox)
      main$getChild()$setShadowType(GtkShadowType["none"])
      
      gSignalConnect(toggleButton , "toggled", .self$onSelectCasesToggled)
      gSignalConnect(button.clear , "clicked", function(...){textBuffer$setText("")})
      gSignalConnect(button.update, "clicked", .self$onSelectCasesUpdated)
    },
    
    onSelectCasesToggled = function(button){
      subset.on <- button$getActive()
      button.update$setSensitive(subset.on)
      button.clear$setSensitive(subset.on)
      textView$setSensitive(subset.on)
      data$setSubset.on(subset.on)
      data$linkDataFrame()
    },
    
    onSelectCasesUpdated = function(button){
      iter <- textBuffer$getBounds()
      text <- textBuffer$getText(iter$start, iter$end)
      text <- localize(text)
      data$setSubset.condition(text)
      result <- data$linkDataFrame()
      if(result) rzTools$getVariableView()$setSubsetSummaries()
    }
    
  ))
selectCases$accessors("main")
