rzPlotScript <- 
setRefClass("RzPlotScript",
  fields = c("script", "freeScript", "aes",
             "theme", "labs", "data", "layercount"),
  methods = list(
    initialize  = function(...) {
      initFields(...)
      script   <<- list()
      freeScript <<- list()
      aes      <<- ""
      labs <<- data.frame(1:2, row.names=c("vars", "labels"), stringsAsFactors=FALSE)
      labs[[1]] <<- NULL
      theme    <<- ""
      layercount <<- numeric()
      data     <<- NULL
    },
    
    setScript = function(layer, type=NULL, args=NULL, add=FALSE){
      if(is.na(layercount[layer])) layercount[layer] <<- 1
      
      if (add) {
        ind <- length(script[[layer]]) + 1
        script[[layer]][[ind]] <<- list(type=type, args=args)
      } else {
        script[[layer]][[layercount[layer]]] <<- list(type=type, args=args)
      }
    },
    
    clearScript = function(layer){
      if(is.na(layercount[layer])) layercount[layer] <<- 1
      if(length(script[[layer]]) < layercount[layer] ) return()
      
      script[[layer]][ seq(layercount[layer],length(script[[layer]])) ] <<- NULL
    },
    
    setFreeScript = function(name, script){
      freeScript[name] <<- script
    },
        
    stackLayer = function(layer){
      if(is.na(layercount[layer])) layercount[layer] <<- 1
      if(length(script[[layer]]) < layercount[layer] ) return()
      
      layercount[layer] <<- length(script[[layer]]) + 1
    },
    
    removeLayer = function(layer){
      if(is.na(layercount[layer])) layercount[layer] <<- 1
      if(layercount[layer] == 1 ) return()
      
      layercount[layer] <<- layercount[layer] - 1
      .self$clearScript(layer)
    },

    getScript = function(){
      text <- character()
      
      data.tmp <- ""
      if (!is.null(data)) {
        data.tmp <- sprintf("data=%s", data$getDataSetName())          
      }
      
      if (any(nzchar(aes))) {
        aes.tmp <- aes[nzchar(aes)]
        aes.tmp.names <- names(aes.tmp)
        aes.tmp <- paste(aes.tmp.names, aes.tmp, sep="=", collapse=", ")
        aes.tmp <- sprintf("aes(%s)", aes.tmp)

        text["data"] <- sprintf("p <- ggplot(%s, %s)", data.tmp, aes.tmp)
      } else {
        text["data"] <- sprintf("p <- ggplot(%s)", data.tmp)
      }
      if (any(nzchar(theme))) {
        theme.tmp <- theme[nzchar(theme)]
        theme.tmp.names <- names(theme.tmp)
        theme.tmp <- paste(theme.tmp.names, theme.tmp, sep="=", collapse=", ")
        text["theme"] <- sprintf("theme(%s)", theme.tmp)
      }

      if (ncol(labs) > 0) {
        labs.tmp <- labs
        title <- labs.tmp[["title"]][2]
        if (is.null(title) || !nzchar(title)) title <- NULL
        labs.tmp[["title"]] <- NULL
        variableNames  <- data$getVariableNames()
        variableLabels <- data$getVariableLabels()
        names(variableLabels) <- variableNames
        
        vals    <- as.character(labs.tmp[1,])
        names(vals) <- labs.tmp.names <- colnames(labs.tmp)
        vals    <- vals[nzchar(vals)]
        labels  <- as.character(labs.tmp[2,])
        names(labels) <- labs.tmp.names
        labels  <- labels[names(vals)]
        labels2 <- variableLabels[vals]
        labs.tmp.names <- names(vals)
        
        labs.tmp <- ifelse(labels  ==gettext("variable name") , vals   , labels)
        labs.tmp <- ifelse(labs.tmp==gettext("variable label"), labels2, labs.tmp)
        names(labs.tmp) <- labs.tmp.names
        if (!is.null(title)) labs.tmp["title"] <- title
        labs.tmp <- labs.tmp[c("title", "x", "y", "group", "fill", "colour", "shape", "size", "linetype", "alpha")]
        labs.tmp <- labs.tmp[!is.na(labs.tmp)]
        labs.tmp.names <- names(labs.tmp)
        
        labs.tmp <- sprintf('%s="%s"', labs.tmp.names, labs.tmp)
        labs.tmp <- paste(labs.tmp, collapse=",")
        if(!nzchar(labs.tmp)) labs.tmp <- NULL
        if(length(labs.tmp) > 0 && nzchar(labs.tmp)) {
          text["labs"] <- sprintf("labs(%s)", labs.tmp)          
        }
      }      
      
      script.tmp <- script
      
      geom         <- .self$buildLayer("geom"        , script.tmp[["geom"]])
      stat         <- .self$buildLayer("stat"        , script.tmp[["stat"]])
      facet        <- .self$buildLayer("facet"       , script.tmp[["facet"]])
      scale_x      <- .self$buildLayer("scale_x"     , script.tmp[["scale_x"]])      
      scale_y      <- .self$buildLayer("scale_y"     , script.tmp[["scale_y"]])      
      coord        <- .self$buildLayer("coord"       , script.tmp[["coord"]])
      scale_fill   <- .self$buildLayer("scale_fill"  , script.tmp[["scale_fill"]])      
      scale_colour <- .self$buildLayer("scale_colour", script.tmp[["scale_colour"]])
      xlim         <- freeScript["xlim"]
      ylim         <- freeScript["ylim"]
      
      text <- c(text["data"], geom, stat, facet, scale_x, scale_y, xlim, ylim, coord,
                scale_fill, scale_colour, text["labs"], "theme_rz()", text["theme"])
      text <- text[!is.na(text)]

      
      text <- paste(text, collapse="\np <- p +")
      text <- c(text, "print(p)")
      text <- tidy.source(text=text, output=FALSE, keep.blank.line=FALSE, width.cutoff=60)$text.tidy
      text <- paste(text, collapse="\n")
      return(text)
    },
    
    setAes = function(type, value){
      aes[type] <<- value
    },
    
    getAes = function(){
      return(aes)
    },
    
    setLabs = function(df){
      labs[names(df)] <<- df
    },

    setTheme = function(type, value){
      theme[type] <<- value
    },
    
    setData = function(data) {
      script <<- list()
      aes      <<- ""
      labs <<- data.frame(1:2, row.names=c("vars", "labels"))
      labs[[1]] <<- NULL
      theme    <<- ""
      data     <<- data
    },
    
    buildLayer = function(layer, list, deparse=TRUE){
      layers <- lapply(list, function(x) {
        args <- x$args
        args <- args[ ! sapply(args, is.null) ]
        args <- unlist(args)
        args.names <- names(args)
        if (is.null(args.names)) args.names <- rep("", length(args))
        
        args <- ifelse(nzchar(args.names), sprintf("%s=%s", args.names, args), args)
        args <- paste(args, collapse=",")
        return( sprintf("%s_%s(%s)", layer, x$type, args) )
      })
      
      if (length(layers) > 0) {
        return(layers)        
      } else {
        return(NULL)
      }
    }
    
  )
)
