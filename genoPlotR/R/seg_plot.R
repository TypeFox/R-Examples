################################
## Plot class and methods
################################
## plot is a list that contains at least three named elements:
## - a function (func)
## - a list of arguments (arguments) to that function
## - an 'xargs' argument designating which variable(s) give x coordinates
## - a  'yargs' argument designating which variable(y) give y coordinates
## - a ylim argument, giving the boundaries on the y axis. Derived
##   automatically from yargs if not set.
seg_plot <- function(func, args=NULL, xargs=c("x", "x0", "x1", "x2", "v"),
                     yargs=c("y", "y0", "y1", "y2", "h"),
                     ylim=NULL){
  as.seg_plot(list(func=func, args=args, xargs=xargs, yargs=yargs, ylim=ylim))
}
as.seg_plot <- function(seg_plot){
  ## Check  that args is a list and that it has func and args
  if (!is.list(seg_plot))
    stop("list argument should be a list object")
  if (!all(c("func", "args") %in% names(seg_plot)))
    stop("list should have elements func, args")
  ## Check that function is a function
  if (!is.function(seg_plot$func))
    stop("func argument should be a function")
  ## Check args
  if (!(is.null(seg_plot$args) || is.list(seg_plot$args)))
    stop("args should be NULL or a list")
  ## Check that by default, the default.units are "native" and not "npc"
  if (is.null(seg_plot$args$default.units))
    seg_plot$args$default.units = "native"
  ## Check xargs and yargs
  if (is.null(seg_plot$xargs)) {
    seg_plot$xargs <- c("x", "x0", "x1", "x2", "v")
  } else if (!is.character(seg_plot$xargs)){
    stop("xargs should be NULL or a character")
  }
  if (is.null(seg_plot$yargs)) {
    seg_plot$yargs <- c("y", "y0", "y1", "y2", "h")
  } else if (!is.character(seg_plot$yargs)){
    stop("yargs should be NULL or a character")
  }
  ## Check ylim
  if (!is.null(seg_plot$ylim) && !is.numeric(seg_plot$ylim))
    stop("ylim should be NULL or numeric")
  ## Check that columns designated by xargs and yargs are numeric or null
  if (!all(sapply(c(seg_plot$xargs, seg_plot$yargs), function(x)
                  is.numeric(seg_plot$args[[x]]) ||
                  is.null(seg_plot$args[[x]]))))
    stop("All arguments designated by xargs & yargs should be numeric or null")
  ## Set class and return
  class(seg_plot) <- c("seg_plot", "list")
  seg_plot
}
is.seg_plot <- function(seg_plot){
  inherits(seg_plot, "seg_plot")
}
## Tries to trim correctly the seg_plot objects...
## Made difficult because of the complexity of the objects
trim.seg_plot <- function(x, xlim=NULL, ...){
  ## Check that we are in the right class
  if (!is.seg_plot(x))
    stop("x should be a seg_plot object.")
  xlim <- as.numeric(xlim)
  ## If xlim is null, return whole object. Else, trim
  if (!is.null(xlim)){
    if (!is.numeric(xlim)) stop("xlim must be numeric")
    if (length(xlim) != 2) stop("xlim must be length 2")
    mainargs <- names(x$args)
    gpargs <- NULL
    ## Define arguments to filter
    if ("gp" %in% mainargs && inherits(x$args$gp, "gpar")){
      mainargs <- mainargs[mainargs != "gp"]
      gpargs <- names(x$args$gp)
    }
    ##std_xargs <- c("x", "x0", "x1", "x2", "v")
    xargs <- mainargs[mainargs %in% x$xargs]
    if (length(xargs) > 0){
      ## Check that all xargs have the same length
      nrows <- length(x$args[[xargs[1]]])
      if (length(xargs) > 1){
        sapply(xargs[2:length(xargs)], function(sp) {
          if(length(x$args[[sp]]) != nrows)
            stop(paste("Argument", sp, "has a different length as argument",
                       xargs[1]))
        })
      }
      ## Get the correct indexes
      idx <- rep(TRUE, nrows)
      for (i in 1:length(xargs)){
        idx <- idx & x$args[[xargs[i]]] >=
          xlim[1] & x$args[[xargs[i]]] <= xlim[2] 
      }
      ## Filter everything...
      for (arg in mainargs){
        if (length(x$args[[arg]]) == nrows)
          x$args[[arg]] <- x$args[[arg]][idx]
      }
      if (!is.null(gpargs)){
        for (arg in gpargs){
          if (length(x$args$gp[[arg]]) == nrows)
            x$args$gp[[arg]] <- x$args$gp[[arg]][idx]
        }
      }
    }
  }
  x
}
## Get a range from a seg_plot object, allowing margin on each side
nice_ylim.seg_plot <- function(seg_plot, margin=0.05){
  if (!is.seg_plot(seg_plot))
    stop("seg_plot should be a seg_plot object")
  if (!is.numeric(margin))
    stop("margin should be numeric")
  if (length(seg_plot$yargs %in% names(seg_plot$args)) < 1)
    stop("No arguments in seg_plot$args matches seg_plot$yargs")
  flock <- unlist(seg_plot$args[seg_plot$yargs])
  rng <- range(flock)
  span <- diff(rng)
  c(rng[1]-span*margin, rng[2]+span*margin)
}
