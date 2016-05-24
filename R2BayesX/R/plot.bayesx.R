plot.bayesx <- function(x, model = NULL, term = NULL, which = 1L, ask = FALSE, ...)
{
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  pc <- FALSE
  which.match <- c("effect", "coef-samples", "var-samples", "intcpt-samples", 
    "hist-resid", "qq-resid", "scatter-resid", "scale-resid", "scale-samples", 
    "max-acf")
  if(!is.character(which)) {
    if(any(which > 10L))
      which <- which[which <= 10L]
    which <- which.match[which]
  } else which <- which.match[pmatch(which, which.match)]
  if(length(which) > length(which.match) || !any(which %in% which.match))
    stop("argument which is specified wrong!")
  x <- get.model(x, model)
  nx <- length(x)
  if("max-acf" %in% which) {
    samples <- NULL
    for(i in 1L:nx) {
      if(!is.null(x[[i]]$effects)) 
        for(k in 1L:length(x[[i]]$effects))
          if(!is.null(x[[i]]$effects[[k]])) {
              if(!"data.frame" %in% class(x[[i]]$effects[[k]])) {
                for(kk in 1L:length(x[[i]]$effects[[k]])) {
                  samples <- cbind(samples, attr(x[[i]]$effects[[k]][[kk]], "sample"))
                  samples <- cbind(samples, attr(x[[i]]$effects[[k]][[kk]], "variance.sample"))
                }
              } else {
                samples <- cbind(samples, attr(x[[i]]$effects[[k]], "sample"))
                samples <- cbind(samples, attr(x[[i]]$effects[[k]], "variance.sample"))
              }
          }
      if(!is.null(x[[i]]$fixed.effects))
        samples <- cbind(samples, attr(x[[i]]$fixed.effects, "sample"))
      samples <- cbind(samples, attr(x[[i]]$variance, "sample"))
    }
  if(!is.null(samples)) {
    args <- list(...)
    args$x <- samples
    args$max.acf <- TRUE
    do.call("plotsamples", args)
  } else warning("nothing to plot!")
  return(invisible(NULL))
  }
  if((!"effect" %in% which) && (!"coef-samples" %in% which) && (!"max-acf" %in% which)
    && (!"var-samples" %in% which) && (!"intcpt-samples" %in% which)) {
    model.names <- names(x)
    pc <- TRUE
    for(i in 1L:nx)
      which.plots(x[[i]], which, ask, model.names[i], nx, ...)
    if(nx > 1L || length(which) > 1L)		
      par(op)
  } else {
    if(is.null(term) && !ask) {
      nt <- 0L
      for(i in 1L:nx)
        nt <- nt + length(x[[i]]$effects)
      if(nt > 1L && !("intcpt-samples" %in% which) && !ask)
        setmfrow(nt)
      if(("intcpt-samples" %in% which) && nx > 1L && !ask)
        setmfrow(nx)
    } else {
      nt <- neffects(x, term)
      if(!ask) {
        if(nt > 1L && !("intcpt-samples" %in% which) && !ask) 
          setmfrow(nt)
        if(("intcpt-samples" %in% which) && nx > 1L && !ask)
          setmfrow(nx)
      }
    }
    args <- list(...)
    for(i in 1L:nx) {
      if("intcpt-samples" %in% which) {
        if(!is.null(attr(x[[i]]$fixed.effects, "sample"))) {
          pc <- TRUE
          args$x <- attr(x[[i]]$fixed.effects, "sample")[,1L]
          args$selected <- "(Intercept)"
          args$var <- FALSE
          if(ask && i == 1L)
            par(ask = TRUE)
          do.call("plotsamples", args)	
        }
      } else {
        if(is.null(term))
          ts <- 1:length(x[[i]]$effects)
        else
          ts <- term
        ne <- names(x[[i]]$effects)
        if(is.null(ne) || !is.character(ts))
          ne <- 1L:length(x[[i]]$effects)
        for(j in 1L:length(ts)) {
          if(ask && j == 1L)
            par(ask = TRUE)
          if(is.character(ts[j])) {
            tmp <- splitme(ts[j])
            tmp <- resplit(tmp[tmp != " "])
            take <- NULL
            for(jj in 1:length(ne)) {
              ## if(!is.na(pmatch(tmp, ne[jj])))
              if(length(grep(tmp, ne[jj], fixed = TRUE)))
                take <- c(take, jj)
            }
          } else take <- match(ts[j], ne)
          if(length(take) > 0L && length(x[[i]]$effects) > 0L && !is.na(take)) {
            for(takeme in take) {
              args$x <- x[[i]]$effects[[takeme]]
              args$diagnostics <- FALSE
              if("coef-samples" %in% which) {
                args$diagnostics <- 1L
                if(length(ts) > 1L && (is.null(args$all.acf) || 
                  !is.null(args$all.acf) && !args$all.acf))
                  par(ask = TRUE)
              }
              if("var-samples" %in% which) 
                args$diagnostics <- 2L
              if(!is.null(args$x)) {
                args$ask <- ask
                pc <- TRUE
                do.call("plot", args)	
              }
            }
          }
        }
      }	
    }		
  }
  if(!pc)
    warning("there is nothing to plot!")

  return(invisible(NULL))
}

