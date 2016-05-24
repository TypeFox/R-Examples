plot.mob <- function(x, terminal_panel = node_bivplot, tnex = NULL, ...) {
  if(is.null(tnex)) {
    xdat <- try(x@data@get("input"), silent = TRUE)
    tnex <- if(inherits(xdat, "try-error")) 2 else 2 * NCOL(xdat)
  }
  plot.BinaryTree(x, terminal_panel = terminal_panel, tnex = tnex, ...)
}

node_scatterplot <- function(mobobj, which = NULL, col = "black", linecol = "red",
  cex = 0.5, pch = NULL, jitter = FALSE, xscale = NULL, yscale = NULL, ylines = 1.5,
  id = TRUE, labels = FALSE)
{
    ## extract dependent variable
    y <- response(mobobj)
    ynam <- names(y)[1]
    y <- y[[1]]
    if(is.factor(y)) y <- as.numeric(y) - 1 
    if(is.Surv(y)) {
      surv <- TRUE
      if(is.null(pch)) pch <- abs(y[,2] - 1) * 18 + 1
      y <- y[,1]
    } else {
      surv <- FALSE
      if(is.null(pch)) pch <- 1
    }
    y <- as.numeric(y)
    pch <- rep(pch, length.out = length(y))

    if(jitter) y <- jitter(y)

    ## extract regressor matrix
    x <- model.matrix(mobobj@tree$model)
    if(is.null(which)) { which <- if(NCOL(x) > 1) 2:NCOL(x) else 1 }
    xnam <- colnames(x)[which]
    x <- x[,which, drop = FALSE]
    k <- NCOL(x)
    
    if(is.null(xscale)) xscale <- apply(x, 2, function(xi) range(xi) + c(-0.1, 0.1) * diff(range(xi)))
        else xscale <- matrix(xscale)
    if(is.null(yscale)) yscale <- range(y) + c(-0.1, 0.1) * diff(range(y))
         
    ## panel function for scatter plots in nodes
    rval <- function(node) {
    
        ## dependent variable setup
	y <- rep.int(y, node$weights)	
	yhat <- fitted(node$model)
	if(!surv) yhat <- rep.int(yhat, node$weights)
        pch <- rep.int(pch, node$weights)

        ## viewport setup
        top_vp <- viewport(layout = grid.layout(nrow = 2*k, ncol = 3,
                           widths = unit(c(ylines, 1, 1), c("lines", "null", "lines")),  
			   heights = unit(rep.int(c(2.5, 1), k) - c(1.5, rep.int(0, 2*k-1)),
			                  rep.int(c("lines", "null"), k))),
                           width = unit(1, "npc"), 
                           height = unit(1, "npc") - unit(2, "lines"),
			   name = paste("node_scatterplot", node$nodeID, sep = ""))
        pushViewport(top_vp)
        grid.rect(gp = gpar(fill = "white", col = 0))

        ## main title
        top <- viewport(layout.pos.col = 2, layout.pos.row = 1)
        pushViewport(top)
	mainlab <- paste(ifelse(id, paste("Node", node$nodeID, "(n = "), ""),
	                 sum(node$weights), ifelse(id, ")", ""), sep = "")
        grid.text(mainlab)
        popViewport()
	
	for(i in 1:k) {
            plot_vpi <- viewport(layout.pos.col = 2, layout.pos.row = 2*i, xscale = xscale[,i],
	        yscale = yscale, name = paste("node_scatterplot", i, node$nodeID, "plot", sep = ""))
            pushViewport(plot_vpi)
	
            ## regressor setup
	    xi <- rep.int(x[,i], node$weights)
	    oi <- order(xi)

            ## scatterplot
	    grid.points(xi, y, gp = gpar(col = col, cex = cex), pch = pch)
            grid.lines(xi[oi], yhat[oi], default.units = "native", gp = gpar(col = linecol))

            grid.xaxis(at = c(ceiling(xscale[1,i]*10), floor(xscale[2,i]*10))/10)
            grid.yaxis(at = c(ceiling(yscale[1]), floor(yscale[2])))

            if(labels) {
                grid.text(xnam[i], x = unit(0.5, "npc"), y = unit(-2, "lines"))
                grid.text(ynam, y = unit(0.5, "npc"), x = unit(-3, "lines"), rot = 90)
	    }

            grid.rect(gp = gpar(fill = "transparent"))
            upViewport()
        }
	upViewport()
    }
	    
    return(rval)
}
class(node_scatterplot) <- "grapcon_generator"

node_bivplot <- function(mobobj, which = NULL, id = TRUE, pop = TRUE,
  pointcol = "black", pointcex = 0.5,
  boxcol = "black", boxwidth = 0.5, boxfill = "lightgray",
  fitmean = TRUE, linecol = "red",
  cdplot = FALSE, fivenum = TRUE, breaks = NULL,
  ylines = NULL, xlab = FALSE, ylab = FALSE, margins = rep(1.5, 4),
  ...)
{
    ## obtain dependent variable
    y <- response(mobobj)
    if(isTRUE(ylab)) ylab <- names(y)
    if(identical(ylab, FALSE)) ylab <- ""
    if(is.null(ylines)) ylines <- ifelse(identical(ylab, ""), 0, 2)
    y <- y[[1]]

    ## obtain explanatory variables
    X <- try(mobobj@data@get("input"), silent = TRUE)
    
    ## if no explanatory variables: behave like plot.BinaryTree
    if(inherits(X, "try-error")) {
      rval <- switch(class(y)[1],
        "Surv" = node_surv(mobobj, id = id, ...),
        "factor" = node_barplot(mobobj, id = id, ...),
        "ordered" = node_barplot(mobobj, id = id, ...),
        node_boxplot(mobobj, ...))
      return(rval)
    }
    
    ## reverse levels for spine/CD plot
    if(is.factor(y)) y <- factor(y, levels = rev(levels(y)))

    ## number of panels needed
    if(is.null(which)) which <- 1:NCOL(X)
    X <- X[,which,drop=FALSE]
    k <- NCOL(X)
    xlab <- if(!identical(xlab, FALSE)) { if(isTRUE(xlab)) colnames(X) else rep(xlab, length.out = k) }
      else rep("", k)
  
    ## set up appropriate panel functions
    if(is.factor(y)) {
      ## CD plots and spine plots
      ## re-use implementation from vcd package
      if(!requireNamespace("vcd", quietly = TRUE)) 
          stop("the `vcd' package is required for CD plots")
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
        grid.xaxis(at = c(ceiling(xscale[1,i]*10), floor(xscale[2,i]*10))/10)
        grid.yaxis(at = c(ceiling(yscale[1]), floor(yscale[2])))
	grid.rect(gp = gpar(fill = "transparent"))
        if(ylab != "") grid.text(ylab, y = unit(0.5, "npc"), x = unit(-2.5, "lines"), rot = 90)
        if(xlab[i] != "") grid.text(xlab[i], x = unit(0.5, "npc"), y = unit(-2, "lines"))                
        if(pop) popViewport() else upViewport()
      }
      cat_fun <- function(x, y, yfit, i, name, ...) {
        xlev <- levels(x)
        pushViewport(plotViewport(margins = margins, name = name,
	  yscale = yscale, xscale = c(0.3, xscale[2,i]+0.7)))

        for(j in seq(along = xlev)) {
	  by <- boxplot(y[x == xlev[j]], plot = FALSE)
          xl <- j - boxwidth/4
	  xr <- j + boxwidth/4

          ## box & whiskers
          grid.lines(unit(c(xl, xr), "native"), 
                     unit(by$stats[1], "native"), gp = gpar(col = boxcol))
          grid.lines(unit(j, "native"), 
                     unit(by$stats[1:2], "native"), gp = gpar(col = boxcol, lty = 2))
          grid.rect(unit(j, "native"), unit(by$stats[2], "native"), 
                    width = unit(boxwidth, "native"), height = unit(diff(by$stats[2:3]), "native"),
                    just = c("center", "bottom"), 
                    gp = gpar(col = boxcol, fill = boxfill))
          grid.rect(unit(j, "native"), unit(by$stats[3], "native"),
                    width = unit(boxwidth, "native"), 
                    height = unit(diff(by$stats[3:4]), "native"),
                    just = c("center", "bottom"), 
                    gp = gpar(col = boxcol, fill = boxfill))
          grid.lines(unit(j, "native"), unit(by$stats[4:5], "native"), 
                     gp = gpar(col = boxcol, lty = 2))
          grid.lines(unit(c(xl, xr), "native"), unit(by$stats[5], "native"), 
                     gp = gpar(col = boxcol))

          ## outlier
          n <- length(by$out)
          if (n > 0) {
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
        grid.xaxis(at = 1:length(xlev), label = xlev)
        grid.yaxis(at = c(ceiling(yscale[1]), floor(yscale[2])))      

        if(ylab != "") grid.text(ylab, y = unit(0.5, "npc"), x = unit(-3, "lines"), rot = 90)
        if(xlab[i] != "") grid.text(xlab[i], x = unit(0.5, "npc"), y = unit(-2, "lines"))                
        if(pop) popViewport() else upViewport()
      }    
    }    

    rval <- function(node) {
    
      ## dependent variable
      y <- rep(y, node$weights)

      ## set up top viewport
      top_vp <- viewport(layout = grid.layout(nrow = k, ncol = 2,
    			 widths = unit(c(ylines, 1), c("lines", "null")), heights = unit(k, "null")),
    			 width = unit(1, "npc"), height = unit(1, "npc") - unit(2, "lines"),
        		 name = paste("node_mob", node$nodeID, sep = ""))
      pushViewport(top_vp)
      grid.rect(gp = gpar(fill = "white", col = 0))

      ## main title
      top <- viewport(layout.pos.col = 2, layout.pos.row = 1)
      pushViewport(top)
      mainlab <- paste(ifelse(id, paste("Node", node$nodeID, "(n = "), ""),
        	       sum(node$weights), ifelse(id, ")", ""), sep = "")
      grid.text(mainlab, y = unit(1, "npc") - unit(0.75, "lines"))
      popViewport()

      for(i in 1:k) {
        ## get x and y 
    	xi <- rep(X[,i], node$weights)
        o <- ORDER(xi)
        yi <- y[o]
        xi <- xi[o]
        yfit <- rep(fitted(mobobj), node$weights)[o]

        ## select panel
    	plot_vpi <- viewport(layout.pos.col = 2, layout.pos.row = i)
    	pushViewport(plot_vpi)

    	## call panel function
    	if(is.factor(xi)) cat_fun(xi, yi, yfit, i, paste("node_mob", node$nodeID, "-", i, sep = ""), ...)
          else num_fun(xi, yi, yfit, i, paste("node_mob", node$nodeID, "-", i, sep = ""), ...)
    	if(pop) popViewport() else upViewport()
      }
      if(pop) popViewport() else upViewport()
    }
    
    return(rval)
}
class(node_bivplot) <- "grapcon_generator"
