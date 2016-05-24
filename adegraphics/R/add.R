setMethod(
  f = "addtext",
  signature = "ADEgORtrellis",
  definition = function(object, xcoord, ycoord, label, which, plot, pos = -1, ...) {
    
    size <- max(length(xcoord), length(ycoord), length(label))
    xcoord <- rep_len(xcoord, length.out = size)
    ycoord <- rep_len(ycoord, length.out = size)
    labels <- rep_len(label, length.out = size)
    
    ## sorting parameters
    sortparameters <- sortparamADEg(...)
    params <- list()
    params$adepar <- list(plabels = list(srt = 0))
    sortparameters <- modifyList(params, sortparameters, keep.null = TRUE)
    params <- sortparameters$adepar$plabels
    
    if(inherits(object, "ADEg")) {
      xlim <- object@g.args$xlim
      ylim <- object@g.args$ylim
    } else {
      xlim <- object$x.limits
      ylim <- object$y.limits
    }
    
    textadded <- xyplot(0 ~ 0, panel = function(x, y, ...) panel.text(xcoord, ycoord, labels, col = params$col, cex = params$cex, alpha = params$alpha, srt = params$srt), plot = FALSE)
    textadded$call <- call("xyplot", 0 ~ 0, xlim = substitute(xlim), ylim = substitute(ylim), xcoord = substitute(xcoord), ycoord = substitute(ycoord), labels = substitute(labels), 
                           col = params$col, cex = params$cex, alpha = params$alpha, srt = params$srt,
                           panel = function(xcoord, ycoord, labels, col, cex, alpha, srt) panel.text(x = xcoord, y = ycoord, labels = labels, col = col, cex = cex, alpha = alpha, srt = srt))
    
    obj <- superpose(object, textadded, plot = FALSE)
    names(obj) <- c("object", "textadded")
    
    if(plot)
      print(obj)
    invisible(obj)
  })


setMethod(
  f = "addtext",
  signature = "ADEgS",
  definition = function(object, xcoord, ycoord, label, which = 1:length(object), plot = TRUE, pos = -1, ...) {
    
    ngraph <- length(object)
    if(max(which) > ngraph)
      stop("Values in 'which' should be lower than the length of object")
    
    ## sorting parameters
    sortparameters <- sortparamADEg(...)
    params <- list()
    params$adepar <- list(plabels = list(alpha = 1, cex = 1, col = "black", srt = 0))
    sortparameters <- modifyList(params, sortparameters, keep.null = TRUE)
    params <- sortparameters$adepar
    params <- lapply(unlist(params, recursive = FALSE), function(X) rep(X, length.out = length(which)))
    
    if(length(which) == 1) { # if only one subgraph is selected, all the labels are displayed on this unique subgraph
      size <- max(length(xcoord), length(ycoord), length(label))
      xcoord <- rep_len(xcoord, length.out = size)
      ycoord <- rep_len(ycoord, length.out = size)
      labels <- rep_len(label, length.out = size)
      
      for (i in which)
        object[[i]] <- addtext(object[[i]], xcoord, ycoord, labels, ..., which = 1, plot = FALSE)
      
    } else { # if several subgraphs are selected, each label is displayed on one subgraph; there is only one label by subgraph
      xcoord <- rep_len(xcoord, length.out = length(which))
      ycoord <- rep_len(ycoord, length.out = length(which))
      labels <- rep_len(label, length.out = length(which))
      for (i in which)
        object[[i]] <- addtext(object[[i]], xcoord[i], ycoord[i], labels[i], which = 1, plot = FALSE,
                               plabels.alpha = params[[1]][i], plabels.cex = params[[2]][i], plabels.col = params[[3]][i], plabels.srt = params[[4]][i])
    }
    
    obj <- object
    if(plot)
      print(obj)
    invisible(obj)
  })