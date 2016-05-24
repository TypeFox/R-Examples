dataHandler <- 
setRefClass("RzDataHandler",
  fields = c("data.handler", "data.collection", "data.set.list", "data.set.list.combo"),
  methods = list(
    initialize            = function(...) {
      initFields(...)
      data.collection <<- data.collection.obj
      data.set.list <<- gtkListStoreNew("character", "character", "character", "character")
      data.set.list.combo <<- gtkComboBoxNewWithModel(data.set.list)
      dsnames        <- data.collection$getDataSetNames()
      original.names <- data.collection$getOriginalNames()
      for( i in seq_len(data.collection$getLength()) ){
        iter <- data.set.list$append()$iter
        data.set.list$set(iter, 0, dsnames[i], 1, sprintf("(%s)", original.names[i]), 2, "", 3, "")
      }
      renderer1 <- gtkCellRendererText()
      renderer2 <- gtkCellRendererText()
      renderer2$setAlignment(1,0.5)
      data.set.list.combo$setFocusOnClick(FALSE)
      data.set.list.combo$packStart(renderer1)
      data.set.list.combo$packStart(renderer2, expand=FALSE)
      data.set.list.combo$addAttribute(renderer1, "text", 0)
      data.set.list.combo$addAttribute(renderer2, "text", 1)
    },
    
    addData  = function(data){
      data.collection$addData(data)
      iter <- data.set.list$append()$iter
      data.set.list$set(iter, 0, data$getData.set.name(), 1, sprintf("( %s )", data$getOriginal.name()), 2, "", 3, "")
      len <- data.set.list$iterNChildren()
      data.set.list.combo$setActive( len-1 )
    },
    
    removeCurrentData = function(){
      data.set.name <- .self$getCurrentDataSetName()
      ind <- data.set.list.combo$getActive()
      data.collection$removeData(data.set.name)
      rm(list=c(data.set.name, paste(data.set.name, ".ds", sep="")),
         envir=.GlobalEnv)
      iter <- data.set.list.combo$getActiveIter()$iter
      data.set.list$remove(iter)
      return(ind)
    },
    
    changeData = function(ind) {
      data.set.list.combo$setActive(ind)
    },
    
    changeDataSetName = function(data.set.name, new.name){
      new.name <- sub("^([[:space:]]+)([^[:space:]]+)([[:space:]]+)$",
                      "\\2", new.name)
      invalid <- grepl("(^$)|(^[0-9]+)|([]\\[\\^$*?|(){}@!\"#$%&'*+,/:;<=>?~[:space:]-])",
                       new.name)
      if (invalid || length(new.name) == 0 ) return(list(result=FALSE, message=gettext("This dataset name is invalid. Please enter a valid dataset name.")))
      invalid <- any(new.name==data.collection$getDataSetNames())
      if (invalid) return(list(result=FALSE, message=gettext("This dataset name is already exists, please enter another dataset name.")))
      data <- data.collection$getData(data.set.name)
      data$setData.set.name(new.name)
      index <- data.set.list.combo$getActive()
      path  <- gtkTreePathNewFromString(index)
      iter  <- data.set.list$getIter(path)$iter
      data.set.list$set(iter, 0, new.name)
      .self$sync()
      return(list(result=TRUE, message=NULL))
    },
    
    getCurrentDataSetName = function(){localize(data.set.list.combo$getActiveText()) },
    
    getCurrentData        = function(){
      data.set.name <- .self$getCurrentDataSetName()
      if(length(data.set.name)==0) return(NULL)
      data <- data.collection$getData(data.set.name)
      return(data)
    },
    
    getData = function(data.set.name){
      data <- data.collection$getData(data.set.name)
      return(data)
    },
    
    sync = function(data.set.name=NULL){
      if (is.null(data.set.name)) {
        data <- .self$getCurrentData()
      } else {
        data <- .self$getData(data.set.name)
      }
      data$linkDataFrame()
    },
    
    syncAll = function(){
      data.collection$syncAll()
    }
    
))
dataHandler$accessors("data.set.list.combo")
