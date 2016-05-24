plotblock <- function(x, residuals = FALSE, range = c(0.3, 0.3), 
  col.residuals = "black", col.lines = "black", c.select = NULL, 
  fill.select = NULL , col.polygons = NULL, data = NULL,
  shift = NULL, trans = NULL, ...)
{
  if(is.null(x))
    return(invisible(NULL))
  if(inherits(x, "formula")) {
    if(is.null(data))
      data <- environment(x)
    else
      if(is.matrix(data))
        data <- as.data.frame(data)
    if(any(grep("+", as.character(x)[2]))) {
      xch <- as.character(x)
      x <- model.frame(as.formula(paste("~", xch[2L])), data = data)
      x <- cbind(model.frame(as.formula(paste("~", xch[3L])), data = data), x)
    } else x <- model.frame(x, data = data)
    if(ncol(x) < 2L)
      stop("formula is specified wrong!")
  }
  is.bayesx <- grepl(".bayesx", class(x))[1L]
  if(is.data.frame(x))
    x <- df2m(x)
  if(!is.list(x) && !is.matrix(x))
    stop("x must be a matrix!")
  if(!is.list(x) && ncol(x) < 2L)
    stop("x must have at least 2 columns!")
  args <- list(...)
  if(is.null(args$xlab))
    args$xlab <- attr(x, "specs")$term
  if(is.null(args$ylab)) {
    if(is.null(attr(x, "specs")$label))
      args$ylab <- paste("f(", args$xlab, ")", sep = "")
    else
      args$ylab <- attr(x, "specs")$label
  }
  if(!is.null(shift))
    shift <- as.numeric(shift[1])
  if(!is.list(x))
    nc <- ncol(x)
  else
    nc <- ncol(x[[1L]])
  if(is.null(c.select)) {
    if(is.bayesx)
      c.select <- c(1L, 2L, 3L, 4L, 6L, 7L) 
    else
      c.select <- 1L:nc
  }
  if(is.null(c.select))
    c.select <- 1L:nc
  if(length(c.select) > nc)
    c.select <- c.select[1L:nc]
  if(is.null(fill.select))
    if(is.bayesx)
      fill.select <- c(0L, 0L, 1L, 2L, 2L, 1L) 
  if(!is.bayesx && length(fill.select) < nc)
    fill.select <- NULL
  xnam <- attr(x, "specs")$term
  if(is.null(xnam))
    xnam <- "x"  
  partial.resids <- NULL
  if(is.null(range)) {
    dow <- 0.3
    up <- 0.3
  } else {
    dow <- range[1L]
    up <- range[2L]	
    if(is.na(dow))
      dow <- 0.3
    if(is.na(up))
      up <- 0.3
  }
  ylim <- NULL
  if(!is.list(x)) {
    if(is.null(e <- attr(x, "partial.resids")))
      residuals <- FALSE
    xu <- unique(x[,1L])
    n <- length(xu)
    effects <- vector("list", n)
    for(i in 1L:n) {
      effects[[i]] <- x[x[,1L] == xu[i], c.select]
      if(!is.matrix(effects[[i]]))
        effects[[i]] <- matrix(effects[[i]], nrow = 1L)
      if(!is.matrix(effects[[i]]))
        effects[[i]] <- matrix(effects[[i]], nrow = 1L)
      if(!is.null(shift)) effects[[i]][, 2L:ncol(x[[i]])] <- effects[[i]][, 2L:ncol(x[[i]])] + shift
      if(!is.null(trans)) {
        if(!is.function(trans)) stop("argument trans must be a function!")
        for(j in 2:ncol(effects[[i]]))
          effects[[i]][, j] <- trans(effects[[i]][, j])
      }
      ylim <- c(ylim, effects[[i]][, 2L:ncol(effects[[i]])])
      colnames(effects[[i]]) <- rep(paste(xnam, xu[i], sep = ""), ncol(effects[[i]]))
      if(residuals) {
        if(length(pres <- e[e[,1L] == xu[i],])) {
          if(is.null(dim(pres)))
            pres <- matrix(pres, nrow = 1)
          if(!is.null(shift))
            pres[, 2L:ncol(pres)] <- pres[, 2L:ncol(pres)] + shift
          attr(effects[[i]], "partial.resids") <- pres
          ylim <- c(ylim, pres[, 2L:ncol(pres)])
        }
      }
    }
    x <- effects
  } else {
    n <- length(x)	
    for(i in 1L:n) {
      if(residuals && !is.null(pres <- attr(x[[i]], "partial.resids"))) {
        pres <- pres[pres[,1L] != 0 & pres[,1L] != -1,]
        if(!is.matrix(pres))
          pres <- matrix(pres, nrow = 1L)
        pres[,1L] <- i
        ylim <- c(ylim, pres[,2L:ncol(pres)])
      }
      if(is.data.frame(x[[i]]))
        x[[i]] <- df2m(x[[i]])
      cn <- colnames(x[[i]])
      cn <- cn[c.select]
      x[[i]] <- x[[i]][,c.select]
      if(!is.matrix(x[[i]]))
        x[[i]] <- matrix(x[[i]], nrow = 1L)
      x[[i]] <- x[[i]][x[[i]][,1L] != 0 & x[[i]][,1L] != -1,]
      if(!is.matrix(x[[i]]))
        x[[i]] <- matrix(x[[i]], nrow = 1L)
      if(nrow(x[[i]]) == 2L && x[[i]][1L, 1L] == -1)
        x[[i]] <- matrix(x[[i]][2L,], nrow = 1L)
      colnames(x[[i]]) <- cn
      if(residuals) {
        if(!is.null(shift)) {
          if(is.matrix(pres))
            pres[, 2L:ncol(pres)] <- pres[, 2L:ncol(pres)] + shift
          else
            pres <- pres + shift
        }
        attr(x[[i]], "partial.resids") <- pres
      }
      if(!is.null(shift)) x[[i]][, 2L:ncol(x[[i]])] <- x[[i]][, 2L:ncol(x[[i]])] + shift
      if(!is.null(trans)) {
        if(!is.function(trans)) stop("argument trans must be a function!")
        for(j in 2:ncol(x[[i]]))
          x[[i]][, j] <- trans(x[[i]][, j])
      }
      ylim <- c(ylim, x[[i]][,2L:ncol(x[[i]])])
    }
  }
  if(is.null(args$xlim))
    args$xlim <- c(0.5, n + 0.5)
  if(is.null(args$ylim))
    args$ylim <- base::range(ylim, na.rm = TRUE)
  if(is.null(args$xlab))
    args$xlab <- xnam	
  if(!is.null(args$add) && args$add)
    par(new = TRUE)
  graphics::plot(args$xlim, args$ylim, type = "n", axes = FALSE, 
    xlab = args$xlab, ylab = args$ylab, main = args$main)
  args <- set.plot2d.specs(ncol(x[[1L]]) - 1L, args, col.lines, is.bayesx)
  xnames <- NULL
  axn <- rep(NA, n)
  args$specs <- args
  args$residuals <- residuals
  args$range <- range
  args$col.residuals <- col.residuals
  args$fill.select <- fill.select
  if(is.null(col.polygons))
    args$col.polygons <- rep(c("grey70", "grey50"), round(nc / 2))
  else
    args$col.polygons <- col.polygons
  args$pb <- TRUE
  for(i in 1L:n) {
    args$x.co <- i
    args$x <- x[[i]]
    if(!is.null(attr(args$x, "partial.resids")))
      attr(args$x, "partial.resids")[,1L] <- i
    do.call(plot2d.default, delete.args(plot2d.default, args))
    axn[i] <- colnames(x[[i]])[1L]
  }
  if(is.null(args$type))
    box()
  else
    if(args$type != "n")
      box()
  if(is.null(args$axes)) {
    axis(2L)
    axis(1L, at = 1L:n, labels = axn)
  } else
    if(args$axes) {
      axis(2L)
      axis(1L, at = 1L:n, labels = axn)
    }

  return(invisible(NULL))
}

