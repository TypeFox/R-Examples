spine <- function(x, ...)
  UseMethod("spine")

spine.formula <- function(formula, data = list(),
  breaks = NULL, ylab_tol = 0.05, off = NULL,
  main = "", xlab = NULL, ylab = NULL, ylim = c(0, 1), margins = c(5.1, 4.1, 4.1, 3.1),
  gp = gpar(), name = "spineplot", newpage = TRUE, pop = TRUE,
  ...)
{
  ## extract x, y from formula
  mf <- model.frame(formula, data = data)
  if(NCOL(mf) != 2) stop("`formula' should specify exactly two variables")
  y <- mf[,1]
  if(!is.factor(y)) stop("dependent variable should be a factor")
  x <- mf[,2]	
  if(is.null(xlab)) xlab <- names(mf)[2]
  if(is.null(ylab)) ylab <- names(mf)[1]
     
  spine(x, y,
    breaks = breaks, ylab_tol = ylab_tol, off = off,
    main = main, xlab = xlab, ylab = ylab, ylim = ylim, margins = margins,
    gp = gp, name = name, newpage = newpage, pop = pop, ...)
}

spine.default <- function(x, y = NULL,
  breaks = NULL, ylab_tol = 0.05, off = NULL,
  main = "", xlab = NULL, ylab = NULL, ylim = c(0, 1), margins = c(5.1, 4.1, 4.1, 3.1),
  gp = gpar(), name = "spineplot", newpage = TRUE, pop = TRUE,
  ...)
{
  ## either supply a 2-way table (i.e., both y and x are categorical)
  ## or two variables (y has to be categorical - x can be categorical or numerical)
  if(missing(y)) {
    if(length(dim(x)) != 2) stop("a 2-way table has to be specified")
    tab <- x
    x.categorical <- TRUE
    if(is.null(xlab)) xlab <- names(dimnames(tab))[1]
    if(is.null(ylab)) ylab <- names(dimnames(tab))[2]
    xnam <- dimnames(tab)[[1]]
    ynam <- dimnames(tab)[[2]]
    ny <- NCOL(tab)
    nx <- NROW(tab)
  } else {
    if(!is.factor(y)) stop("dependent variable should be a factor")
    x.categorical <- is.factor(x)
    if(!x.categorical) stopifnot(is.numeric(x), is.vector(x))
    if(is.null(xlab)) xlab <- deparse(substitute(x))
    if(is.null(ylab)) ylab <- deparse(substitute(y))
    if(x.categorical) {
      tab <- table(x, y)
      xnam <- levels(x)
      nx <- NROW(tab)
    }
    ynam <- levels(y)
    ny <- length(ynam)
  }

  ## graphical parameters
  if(is.null(gp$fill)) gp$fill <- gray.colors(ny)
  gp$fill <- rep(gp$fill, length.out = ny)
  off <- if(!x.categorical) 0 else if(is.null(off)) 0.02 else off/100

  if(x.categorical) {
    ## compute rectangle positions on x axis
    xat <- c(0, cumsum(prop.table(margin.table(tab, 1)) + off))
  } else {
    ## compute breaks for x
    if(is.null(breaks)) breaks <- list()
    if(!is.list(breaks)) breaks <- list(breaks = breaks)
    breaks <- c(list(x = x), breaks)
    breaks$plot <- FALSE
    breaks <- do.call("hist", breaks)$breaks
    ## categorize x
    x1 <- cut(x, breaks = breaks, include.lowest = TRUE)
    ## compute rectangle positions on x axis
    xat <- c(0, cumsum(prop.table(table(x1))))
    ## construct table
    tab <- table(x1, y)
    nx <- NROW(tab)
  }
  
  ## compute rectangle positions on y axis
  yat <- rbind(0, apply(prop.table(tab, 1), 1, cumsum))
  
  ## setup plot
  if(newpage) grid.newpage()
  pushViewport(plotViewport(xscale = c(0, 1 + off * (nx-1)), yscale = ylim,
    default.units = "native", name = name, margins = margins, ...))

  ## compute coordinates
  ybottom <- as.vector(yat[-(ny+1),])
  ybottom[ybottom < ylim[1]] <- ylim[1]
  ybottom[ybottom > ylim[2]] <- ylim[2]
  ytop <- as.vector(yat[-1,])
  ytop[ytop < ylim[1]] <- ylim[1]
  ytop[ytop > ylim[2]] <- ylim[2]
  xleft <- rep(xat[1:nx], rep(ny, nx))
  xright <- rep(xat[2:(nx+1)] - off, rep(ny, nx))
  gp$fill <- rep(gp$fill, nx)

  ## plot rectangles
  grid.rect(xleft, ybottom, width = (xright-xleft), height = (ytop-ybottom),
    just = c("left", "bottom"), default.units = "native", gp = gp)

  ## axes
  ## 1: either numeric or level names
  if(x.categorical)
    grid.text(x = unit((xat[1:nx] + xat[2:(nx+1)] - off)/2, "native"), y = unit(-1.5, "lines"),
      label = xnam, check.overlap = TRUE)
  else
    grid.xaxis(at = xat, label = breaks)
    
  ## 2: axis with level names of y
  yat <- yat[,1]
  equidist <- any(diff(yat) < ylab_tol)
  yat <- if(equidist) seq(1/(2*ny), 1-1/(2*ny), by = 1/ny)
           else (yat[-1] + yat[-length(yat)])/2
  grid.text(x = unit(-1.5, "lines"), y = unit(yat, "native"), label = ynam,
    rot = 90, check.overlap = TRUE)
  
  ## 3: none
  ## 4: simple numeric  
  grid.yaxis(main = FALSE)
  
  ## annotation
  grid.text(xlab, y = unit(-3.5, "lines"))
  grid.text(ylab, x = unit(-3, "lines"), rot = 90)
  grid.text(main, y = unit(1, "npc") + unit(2, "lines"), gp = gpar(fontface = "bold"))

  ## pop
  if(pop) popViewport()
  
  ## return table visualized
  names(dimnames(tab)) <- c(xlab, ylab)
  invisible(tab)  
}
