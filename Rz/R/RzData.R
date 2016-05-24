rzdata <-
setRefClass("RzData",
  fields = c("version", "file.path", "data.set.name", "original.name",
             "data.set", "original.data.set", "data.frame",
             "data.set.subset", "data.frame.subset",
             "log", "subset.condition", "subset.on",
             "system", "package.version", "encoding", "time", "autosaved"),
  methods = list(
    initialize        = function(...) {
      initFields(...)
      version <<- numeric_version("1")
#      data.set.name    <<- NULL
      original.data.set <<- data.set
      data.frame        <<- suppressWarnings(as.data.frame(data.set))
#      log              <<- NULL
      system            <<- R.Version()
      package.version   <<- list(Rz = packageVersion("Rz"),
                                 memisc=packageVersion("memisc"),
                                 RGtk2=packageVersion("RGtk2"))
      encoding          <<- localeToCharset()
      time              <<- Sys.time()
      subset.condition  <<- ""
      subset.on         <<- FALSE
      data.set.subset   <<- NULL
      data.frame.subset <<- NULL
      autosaved         <<- FALSE
    },

    save = function(file){
      rzdata <- new.env()
      rzdata$file.path           <- file.path
      rzdata$data.set.name       <- data.set.name
      rzdata$original.name       <- original.name
      rzdata$data.set            <- data.set
#      rzdata$original.data.set   <- original.data.set
#      rzdata$data.frame          <- data.frame
      rzdata$log                 <- log
      rzdata$system              <- system
      rzdata$package.version     <- package.version
      rzdata$encoding            <- encoding
      rzdata$time                <- time
      base::save(rzdata, file=file)
    },

    revert = function(){
      autosaved  <<- FALSE
      data.set   <<- original.data.set
      data.frame <<- suppressWarnings(as.data.frame(data.set))
      .self$linkDataFrame()
    },

    reloadFromGlobalEnv = function(){
      autosaved  <<- FALSE
      data.set.tmp <- try(get(paste(data.set.name, ".ds", sep=""), envir=.GlobalEnv), silent=TRUE)
      if(is.data.set(data.set.tmp)){
        data.set   <<- data.set.tmp
        data.frame <<- suppressWarnings(as.data.frame(data.set))
        .self$linkDataFrame()
        return(TRUE)
      } else {
        return(data.set.name)
      }
    },
    
    existVar = function(varname){
      var <- data.set[[varname]]
      return(!is.null(var))
    },
    
    addItem = function(item, name){
      autosaved  <<- FALSE
      data.set[[name]]  <<- item
      data.frame        <<- suppressWarnings(as.data.frame(data.set))      
      .self$linkDataFrame()
    },
    
    deleteVars = function(inds){
      autosaved  <<- FALSE
      inds <- sort(inds)
      inds <- rev(inds)
      for(i in inds){
        data.set[[i]] <<- NULL
      }
      data.frame <<- suppressWarnings(as.data.frame(data.set))
      .self$linkDataFrame()      
    },
    
    duplicate = function(inds){
      autosaved  <<- FALSE
      dup <- data.set[inds]
      orig.names <- colnames(data.set)
      dup.names  <- colnames(dup)
      dup.names  <- dup.names.tmp <- sprintf("%s.copy", dup.names)
      for(i in seq_along(dup.names)){
        suffix <- 2
        while(any(dup.names[i]==orig.names)){
          dup.names[i] <- sprintf("%s%s", dup.names.tmp[i], suffix)
          suffix <- suffix + 1
        }
      }
      data.set <<- cbind(data.set, dup)
      colnames(data.set) <<- c(orig.names, dup.names)
      data.frame <<- suppressWarnings(as.data.frame(data.set))
      .self$linkDataFrame()      
    },
    
    ncol              = function() { length(description(data.set)) },
    
    constructVariable = function(col){
      data.frame[[col]] <<- suppressWarnings(as.data.frame(data.set[[col]])[[1]])
      names(data.frame) <<- names(data.set)
    },
    
    linkDataFrame     = function(){
      if(subset.on & nzchar(subset.condition)){
        data.set.subset <<- tryCatch(eval(parse(text=paste("subset(data.set, subset=", subset.condition, ")"))), error=function(...) return(FALSE))
        if(!is.data.set(data.set.subset)){
          dialog <- gtkMessageDialogNew(rzTools$getWindow(), "destroy-with-parent",
                                        "error", "close", gettext("Condition of selecting cases is invalid."))
          dialog$run()
          dialog$destroy()
          return(FALSE)
        }
        data.frame.subset <<- eval(parse(text=paste("subset(data.frame, subset=", subset.condition, ")")))
        rzTools$sync(paste(data.set.name, ".ds", sep=""), data.set.subset)
        rzTools$sync(data.set.name                      , data.frame.subset)
      } else {
        rzTools$sync(paste(data.set.name, ".ds", sep=""), data.set)
        rzTools$sync(data.set.name                      , data.frame)
      }
      if (rzSettings$getAutosave()) saveSession()
      return(TRUE)
    },

    import = function(data){
      
    },

    # getter methods
    getVariableLabels = function() gsub("(^')|('$)", "", as.character(description(data.set))),
    
    getVariableNames  = function() names(description(data.set)),
    
    getMeasurement    = function() sapply(data.set, measurement),
    
    getMissingValues  = function() {
      miss.val  <- lapply(data.set, missing.values)
      miss.val  <- sapply(miss.val, function(x){
        if      (is.null(x))        return("")
        else if (is.null(x@filter)) return(sprintf("range=c(%s, %s)", x@range[1], x@range[2]))
        else                        return(paste(x@filter, collapse=","))
      })
    },
    
    getValueLabels    = function() {
      val.labs  <- lapply(data.set, labels)
      val.labs  <- sapply(val.labs, function(x){ifelse(is.null(x), "", paste(x@values, " \"",x@.Data, "\"", sep="", collapse=", "))})
    },
    
    getSummaries      = function(subset=FALSE) {
      summaries <- NULL
      if (subset) {
        summaries <- lapply(data.set.subset, summary)
      } else {
        summaries <- lapply(data.set, summary)        
      }
      val <- lapply(summaries, as.character)
      lab <- lapply(summaries, names)
      summaies <- sapply(seq_along(val),
                           function(x){
                            lnchar <- nchar(sub("", "", lab[[x]]), type="width") # "sub()" to remove invalid character
                            lnchar <- max(lnchar) - lnchar
                            vnchar <- nchar(val[[x]], type="width")
                            vnchar <- max(vnchar) - vnchar
                            paste(lab[[x]],
                                  sapply(sapply(lnchar, function(x) rep(" ", x)), paste, collapse=""),
                                  "\t:",
                                  sapply(sapply(vnchar, function(x) rep(" ", x)), paste, collapse=""),
                                  val[[x]], sep="", collapse="\n")
                           }
                         )
    },
    
    getSummary = function(row, subset=FALSE){
      summary <- NULL
      if (subset) {
        summary <- summary(data.set.subset[[row]])        
      } else {
        summary <- summary(data.set[[row]])
      }
      lab <- names(summary)
      val <- as.character(summary)
      lnchar <- nchar(sub("", "", lab), type="width") # "sub()" is for removing invalid character.
      lnchar <- max(lnchar) - lnchar
      vnchar <- nchar(val, type="width")
      vnchar <- max(vnchar) - vnchar
      summary <- paste(lab,
                       sapply(sapply(lnchar, function(x) rep(" ", x)), paste, collapse=""),
                       "\t:",
                       sapply(sapply(vnchar, function(x) rep(" ", x)), paste, collapse=""),
                       val, sep="", collapse="\n")
      return(summary)
    },
    
    getDataSetName = function(){
      return(data.set.name)
    }
  )
)
rzdata$accessors(c("version", "data.set.name", "original.name", "data.set", "data.frame", "file.path",
                   "subset.condition", "subset.on", "data.set.subset", "data.frame.subset", "autosaved"))
