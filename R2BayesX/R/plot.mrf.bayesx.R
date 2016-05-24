plot.mrf.bayesx <- function(x, diagnostics = FALSE, ...)
{
  if(!is.null(x)) {
		args <- list(...)
    if(diagnostics == FALSE) {
      if(inherits(x, "random.bayesx")) {
        cx <- colnames(x)
        if(!is.null(args$total) && args$total && any(grep("_tot", cx))) {
          xattr <- attributes(x)
          id <- c(1L, grep("_tot", cx))
          x <- x[,id]
          for(na in names(xattr))
            if(na != "dim" && na != "dimnames" && na != "names")
              attr(x, na) <- eval(parse(text = paste("xattr$", na, sep = "")))
          args$total <- NULL
        } else {
          xattr <- attributes(x)
          id <- 1L:ncol(x)
          id <- id[!grepl("tot_sim", cx) & !grepl("_tot", cx) & !grepl("_sim", cx)]
          x <- x[,id]
          for(na in names(xattr))
            if(na != "dim" && na != "dimnames" && na != "names")
              attr(x, na) <- eval(parse(text = paste("xattr$", na, sep = "")))
        }
      }
      if(is.logical(args$map))
        domap <- args$map
      else {
        if(is.null(args$map))
          domap <- FALSE
        else
          domap <- TRUE
      }
      if(!is.null(attr(x, "map.name")) && domap) {        
        tmp <- try(eval(parse(text = attr(x, "map.name")), envir = globalenv()), silent = TRUE)
        if(all(class(tmp) %in% c("bnd", "list")))
          args$map <- tmp
        if(all(class(args$map) == "try-error"))
          args$map <- NULL
      }
      domap <- FALSE
      if(!is.null(args$map)) {
        if(is.logical(args$map))
          domap <- args$map
        else domap <- TRUE
      } else domap <- FALSE
      if(!domap) {
        args$x <- x
        args$x <- do.call(compute.x.id, delete.args(compute.x.id, args))$x
        if(is.null(args$from))
          args$from <- min(args$x, na.rm = TRUE)
        if(is.null(args$to))
          args$to <- max(args$x, na.rm = TRUE)
        if(is.null(args$main))
          args$main <- paste("Kernel density estimate of term", attr(x, "specs")$label)
        args2 <- args
        args2$x <- do.call(stats::density, delete.args(stats::density.default, args))
        do.call(plot, delete.args("plot.density", args2, package = "stats"))
        if(!is.null(args$kde.quantiles))
          if(args$kde.quantiles) {
            if(is.null(args$lty))
              args$lty <- 2L
            if(is.null(args$col))
              args$col <- "black"
            abline(v = do.call(kde.quantiles, 
              delete.args(kde.quantiles, args)), 
              lty = args$lty, col = args$col)
          }       
      } else {
        if(is.logical(args$map))
          stop("could not find the corresponding map!")
        args$x <- x
        do.call("plotmap", args)
      }
    } else coeffun(x, args, diagnostics)
  } else warning("there is nothing to plot!")

  return(invisible(NULL))
}

