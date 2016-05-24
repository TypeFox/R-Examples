data.collection <-
setRefClass("RzDataCollection",
  fields = c("data.collection"),
  methods = list(
    initialize       = function(...) {
      data.collection <<- list()
      initFields(...)
    },
    
    saveSession = function(file=NULL, force=FALSE){
      if (is.null(file)) {
        autosaved <- sapply(data.collection, function(x) x$getAutosaved())
        if(!force && all(autosaved)) return()
        sapply(data.collection, function(x) x$setAutosaved(TRUE))
        if(!checkConfDir()) return()
        file <- file.path(rzConfPath(), "session.rzs")
      }
      
      session <- data.collection
      save(session, file=file)
    },

    loadSession = function(file=NULL){
      if (is.null(file)) file <- file.path(rzConfPath(), "session.rzs")
      if (file.exists(file)) {
        load(file=file)
        data.collection <<- session
      }
    },
    
    syncAll = function(){
      lapply(data.collection, function(x) x$linkDataFrame())
    },

    addData          = function(data) {
      names <- .self$getDataSetNames()
      match <- grep("^dataset*", names, value=TRUE)
      match <- suppressWarnings(as.numeric(gsub("dataset", "", match)))
      match <- c(0, match)
      data$setData.set.name(sprintf("dataset%s", max(match, na.rm=TRUE)+1 ))
      data.collection[[ length(data.collection) + 1 ]] <<- data
    },

    removeData       = function(data.set.name){
      names <- .self$getDataSetNames()      
      data.collection[ which(names == data.set.name) ] <<- NULL
      if (rzSettings$getAutosave()) {
        saveSession(force=TRUE)
      }
    },

    getDataSetNames  = function(){
      names <- lapply(data.collection, function(x) x$getData.set.name())
      return( unlist(names) )
    },

    getOriginalNames = function(){
      names <- lapply(data.collection, function(x) x$getOriginal.name())
      return( unlist(names) )
    },

    getData          = function(data.set.name){
      names <- .self$getDataSetNames()      
      data  <- data.collection[[ which(names == data.set.name) ]]
      return(data)
    },

    getLength        = function(){length(data.collection)}
  )
)
data.collection$accessors(c("data.collection"))
