completion.combo <- 
setRefClass("RzCompletionCombo",
              fields = c("combo", "entry.completion", "width", "show.combo"),
  methods = list(
    initialize  = function(...) {
      initFields(...)
      if (class(width) == "uninitializedField")      width      <<- 120
      if (class(show.combo) == "uninitializedField") show.combo <<- TRUE
      
      entry.completion <<- gtkEntryCompletionNew()
      entry.completion$setTextColumn(column.definition["vars"])
      entry.completion$setInlineCompletion(TRUE)
      entry.completion$setInlineSelection(TRUE)
      entry.completion$setPopupSetWidth(FALSE)
      entry.completion$clear()
      renderer1 <- gtkCellRendererPixbuf()
      renderer2 <- gtkCellRendererText()
      renderer3 <- gtkCellRendererText()
      entry.completion$packStart(renderer1, expand=FALSE)
      entry.completion$packStart(renderer2, expand=FALSE)
      entry.completion$packStart(renderer3, expand=TRUE)
      entry.completion$addAttribute(renderer1, "pixbuf", column.definition["msr.image"])
      entry.completion$addAttribute(renderer2, "text"  , column.definition["vars"])
      entry.completion$addAttribute(renderer3, "text"  , column.definition["var.labs"])
      
      combo <<- gtkComboBoxEntryNew(show=show.combo)
      combo["width-request"] <<- width
      entry <- combo$getChild()
      entry$setCompletion(entry.completion)
      combo["text-column"] <<- column.definition["vars"]
      combo$clear()
      renderer1 <- gtkCellRendererPixbuf()
      renderer2 <- gtkCellRendererText()
      renderer3 <- gtkCellRendererText()
      combo$packStart(renderer1, expand=FALSE)
      combo$packStart(renderer2, expand=FALSE)
      combo$packStart(renderer3, expand=TRUE)
      combo$addAttribute(renderer1, "pixbuf", column.definition["msr.image"])
      combo$addAttribute(renderer2, "text"  , column.definition["vars"])
      combo$addAttribute(renderer3, "text"  , column.definition["var.labs"])
    },
    
    setText = function(txt){
      entry <- combo$getChild()
      entry$setText(txt)
    },
    
    clear = function(){
      entry <- combo$getChild()
      entry$setText("")
      entry$setText("")
    },
    
    getActiveText = function(){
      txt <- combo$getActiveText()
      return(txt)
    },
    
    setModel = function(model){
      entry.completion$setModel(model)
      if(!is.null(model)){
        model.filterd <- gtkTreeModelFilterNew(model)
        model.filterd$setVisibleColumn(column.definition["select"])
        combo$setModel(model.filterd)
      }
    }
  )
)
completion.combo$accessors("combo")
