node_bivplot <- function(mobobj, which = NULL, id = TRUE, pop = TRUE,
  pointcol = "black", pointcex = 0.5,
  boxcol = "black", boxwidth = 0.5, boxfill = "lightgray", bg = "white",
  fitmean = TRUE, linecol = "red",
  cdplot = FALSE, fivenum = TRUE, breaks = NULL,
  ylines = NULL, xlab = FALSE, ylab = FALSE, margins = rep(1.5, 4),
  mainlab = NULL, ...)
{
  ## obtain dependent variable
  mf <- model.frame(mobobj)
  y <- Formula::model.part(mobobj$info$Formula, mf, lhs = 1L, rhs = 0L)
  if(isTRUE(ylab)) ylab <- names(y)
  if(identical(ylab, FALSE)) ylab <- ""
  if(is.null(ylines)) ylines <- ifelse(identical(ylab, ""), 0, 2)
  y <- y[[1L]]

  ## obtain explanatory variables
  X <- Formula::model.part(mobobj$info$Formula, mf, lhs = 0L, rhs = 1L)
  
  ## fitted node ids
  fitted <- mobobj$fitted[["(fitted)"]]
  
  ## if no explanatory variables: behave like plot.constparty
  if(inherits(X, "try-error")) {
    rval <- switch(class(y)[1L],
      "Surv" = node_surv(mobobj, id = id, mainlab = mainlab, ...),
      "factor" = node_barplot(mobobj, id = id, mainlab = mainlab, ...),
      "ordered" = node_barplot(mobobj, id = id, mainlab = mainlab, ...),
      node_boxplot(mobobj, ...))
    return(rval)
  }
  
  ## reverse levels for spine/CD plot
  if(is.factor(y)) y <- factor(y, levels = rev(levels(y)))

  ## number of panels needed
  if(is.null(which)) which <- 1L:NCOL(X)
  X <- X[,which,drop=FALSE]
  k <- NCOL(X)
  xlab <- if(!identical(xlab, FALSE)) { if(isTRUE(xlab)) colnames(X) else rep(xlab, length.out = k) }
    else rep("", k)

  ## set up appropriate panel functions
  if(is.factor(y)) {
    ## CD plots and spine plots
    ## re-use implementation from vcd package
    if(!requireNamespace("vcd")) stop(sprintf("Package %s is required for spine/CD plots", sQuote("vcd")))
    if(cdplot) {
      num_fun <- function(x, y, yfit, i, name, ...) {
	vcd::cd_plot(x, y, xlab = xlab[i], ylab = ylab, name = name, newpage = FALSE,
          margins = margins, pop = FALSE, ...)
	if(fitmean) {
          #FIXME# downViewport(name = name)
	  grid.lines(x, yfit, default.units = "native", gp = gpar(col = linecol))
          if(pop) popViewport() else upViewport()
        } else {
          #FIXME#
          if(pop) popViewport() else upViewport()
        }
      }
    } else {	   
      xscale <- if(is.null(breaks)) {
        if(fivenum) lapply(X, function(z) {if(is.factor(z)) 1 else fivenum(z) })
	  else lapply(X, function(z) {if(is.factor(z)) 1 else hist(z, plot = FALSE)$breaks })
      } else {
        if(is.list(breaks)) breaks else list(breaks)
      }
      num_fun <- function(x, y, yfit, i, name, ...) {
	vcd::spine(x, y, xlab = xlab[i], ylab = ylab, name = name, newpage = FALSE,
          margins = margins, pop = FALSE, breaks = xscale[[i]], ...)
	if(fitmean) {
          #FIXME# downViewport(name = name)
          xaux <- cut(x, breaks = xscale[[i]], include.lowest = TRUE)	  
          yfit <- unlist(tapply(yfit, xaux, mean))
          xaux <- prop.table(table(xaux))
          xaux <- cumsum(xaux) - xaux/2
	  grid.lines(xaux, yfit, default.units = "native", gp = gpar(col = linecol))
	  grid.points(xaux, yfit, default.units = "native",
            gp = gpar(col = linecol, cex = pointcex), pch = 19)
          if(pop) popViewport() else upViewport()
        } else {
          #FIXME#
          if(pop) popViewport() else upViewport()
        }
      }
    }
    cat_fun <- function(x, y, yfit, i, name, ...) {
      vcd::spine(x, y, xlab = xlab[i], ylab = ylab, name = name, newpage = FALSE,
        margins = margins, pop = FALSE, ...)
      if(fitmean) {
        #FIXME# downViewport(name = name)
        yfit <- unlist(tapply(yfit, x, mean))
        xaux <- prop.table(table(x))
        xaux <- cumsum(xaux + 0.02) - xaux/2 - 0.02
	grid.lines(xaux, yfit, default.units = "native", gp = gpar(col = linecol))
	grid.points(xaux, yfit, default.units = "native",
          gp = gpar(col = linecol, cex = pointcex), pch = 19)
        if(pop) popViewport() else upViewport()
      } else {
        #FIXME#
        if(pop) popViewport() else upViewport()
      }
    }
  } else {
    xscale <- sapply(X, function(z) {if(is.factor(z)) c(1, length(levels(z))) else range(z) })
    yscale <- range(y) + c(-0.1, 0.1) * diff(range(y))

    ## scatter plots and box plots
    num_fun <- function(x, y, yfit, i, name, ...) {
      xscale[,i] <- xscale[,i] + c(-0.1, 0.1) * diff(xscale[,i])
      pushViewport(plotViewport(margins = margins, name = name,
        yscale = yscale, xscale = xscale[,i]))
      grid.points(x, y, gp = gpar(col = pointcol, cex = pointcex))
      if(fitmean) {   
	grid.lines(x, yfit, default.units = "native", gp = gpar(col = linecol))
      }
      grid.xaxis(at = c(ceiling(xscale[1L,i]*10), floor(xscale[2L,i]*10))/10)
      grid.yaxis(at = c(ceiling(yscale[1L]), floor(yscale[2L])))
      grid.rect(gp = gpar(fill = "transparent"))
      if(ylab != "") grid.text(ylab, y = unit(0.5, "npc"), x = unit(-2.5, "lines"), rot = 90)
      if(xlab[i] != "") grid.text(xlab[i], x = unit(0.5, "npc"), y = unit(-2, "lines")) 	       
      if(pop) popViewport() else upViewport()
    }
    cat_fun <- function(x, y, yfit, i, name, ...) {
      xlev <- levels(x)
      pushViewport(plotViewport(margins = margins, name = name,
        yscale = yscale, xscale = c(0.3, xscale[2L,i]+0.7)))

      for(j in seq(along = xlev)) {
        by <- boxplot(y[x == xlev[j]], plot = FALSE)
	xl <- j - boxwidth/4
        xr <- j + boxwidth/4

	## box & whiskers
	grid.lines(unit(c(xl, xr), "native"), 
		   unit(by$stats[1L], "native"), gp = gpar(col = boxcol))
	grid.lines(unit(j, "native"), 
		   unit(by$stats[1L:2L], "native"), gp = gpar(col = boxcol, lty = 2))
	grid.rect(unit(j, "native"), unit(by$stats[2L], "native"), 
		  width = unit(boxwidth, "native"), height = unit(diff(by$stats[2:3]), "native"),
		  just = c("center", "bottom"), 
		  gp = gpar(col = boxcol, fill = boxfill))
	grid.rect(unit(j, "native"), unit(by$stats[3L], "native"),
		  width = unit(boxwidth, "native"), 
		  height = unit(diff(by$stats[3L:4L]), "native"),
		  just = c("center", "bottom"), 
		  gp = gpar(col = boxcol, fill = boxfill))
	grid.lines(unit(j, "native"), unit(by$stats[4L:5L], "native"), 
		   gp = gpar(col = boxcol, lty = 2))
	grid.lines(unit(c(xl, xr), "native"), unit(by$stats[5L], "native"), 
		   gp = gpar(col = boxcol))

	## outlier
	n <- length(by$out)
	if (n > 0L) {
	  grid.points(unit(rep.int(j, n), "native"),  unit(by$out, "native"),
		      size = unit(0.5, "char"), gp = gpar(col = boxcol))
	}	
      }
      if(fitmean) {
        yfit <- unlist(tapply(yfit, x, mean))
	grid.lines(seq(along = xlev), yfit, default.units = "native", gp = gpar(col = linecol))
	grid.points(seq(along = xlev), yfit, default.units = "native",
          gp = gpar(col = linecol, cex = pointcex), pch = 19)
      }
      grid.rect(gp = gpar(fill = "transparent"))
      grid.xaxis(at = 1L:length(xlev), label = xlev)
      grid.yaxis(at = c(ceiling(yscale[1L]), floor(yscale[2L])))      

      if(ylab != "") grid.text(ylab, y = unit(0.5, "npc"), x = unit(-3, "lines"), rot = 90)
      if(xlab[i] != "") grid.text(xlab[i], x = unit(0.5, "npc"), y = unit(-2, "lines")) 	       
      if(pop) popViewport() else upViewport()
    }	 
  }    

  rval <- function(node) {
    
    ## node index
    nid <- id_node(node)
    ix <- fitted %in% nodeids(mobobj, from = nid, terminal = TRUE)

    ## dependent variable
    y <- y[ix]

    ## set up top viewport
    top_vp <- viewport(layout = grid.layout(nrow = k, ncol = 2,
		       widths = unit(c(ylines, 1), c("lines", "null")), heights = unit(k, "null")),
		       width = unit(1, "npc"), height = unit(1, "npc") - unit(2, "lines"),
		       name = paste("node_mob", nid, sep = ""))
    pushViewport(top_vp)
    grid.rect(gp = gpar(fill = bg, col = 0))

    ## main title
    top <- viewport(layout.pos.col = 2, layout.pos.row = 1)
    pushViewport(top)

    if (is.null(mainlab)) { 
      mainlab <- if(id) {
    	function(id, nobs) sprintf("Node %s (n = %s)", id, nobs)
      } else {
    	function(id, nobs) sprintf("n = %s", nobs)
      }
    }
    if (is.function(mainlab)) {
      mainlab <- mainlab(nid, info_node(node)$nobs)
    }
    grid.text(mainlab, y = unit(1, "npc") - unit(0.75, "lines"))
    popViewport()

    for(i in 1L:k) {
      ## get x and y 
      xi <- X[ix, i]
      o <- order(xi)
      yi <- y[o]
      xi <- xi[o]
      yfit <- if(is.null(node$info$object)) {
        fitted(refit.modelparty(mobobj, node = nid))[o]
      } else {
        fitted(node$info$object)[o]
      }

      ## select panel
      plot_vpi <- viewport(layout.pos.col = 2L, layout.pos.row = i)
      pushViewport(plot_vpi)

      ## call panel function
      if(is.factor(xi)) cat_fun(xi, yi, yfit, i, paste("node_mob", nid, "-", i, sep = ""), ...)
	else num_fun(xi, yi, yfit, i, paste("node_mob", nid, "-", i, sep = ""), ...)
      if(pop) popViewport() else upViewport()
    }
    if(pop) popViewport() else upViewport()
  }
  
  return(rval)
}
class(node_bivplot) <- "grapcon_generator"
