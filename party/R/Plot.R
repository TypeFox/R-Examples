
# $Id: Plot.R 512 2013-01-09 10:53:30Z thothorn $

## utility functions for querying the number of
## terminal nodes and the maximal depth of (sub-)trees
nterminal <- function(node) {
    if (node$terminal) return(1)
    nl <- nterminal(node$left)
    nr <- nterminal(node$right)
    return(nl + nr)
}

maxdepth <- function(node) {
    if (node$terminal) return(1)
    nl <- maxdepth(node$left)
    nr <- maxdepth(node$right)
    return(max(c(nl, nr)) + 1)
}


## panel functions for labeling trees:
## inner and terminal nodes and edges.

node_inner <- function(ctreeobj,
                       digits = 3,
		       abbreviate = FALSE,
		       fill = "white",
		       pval = TRUE,
		       id = TRUE)
{
    getLabel1 <- function(x) {
        if (x$terminal) return(rep.int("", 2))
        varlab <- ifelse(abbreviate > 0,
            abbreviate(x$psplit$variableName, as.numeric(abbreviate)),
	    x$psplit$variableName)
	if(pval) {
            pvalue <- 1 - x$criterion$maxcriterion
            plab <- ifelse(pvalue < 10^(-digits),
                           paste("p <", 10^(-digits)),
                           paste("p =", round(pvalue, digits = digits)))
	} else {
	    plab <- ""
	}
        return(c(varlab, plab))
    }

    maxstr <- function(node) {
        lab <- getLabel1(node)
        msl <- ifelse(node$terminal, "", maxstr(node$left))
        msr <- ifelse(node$terminal, "", maxstr(node$right))
        lab <- c(lab, msl, msr)
        return(lab[which.max(nchar(lab))])
    }

    nstr <- maxstr(ctreeobj@tree) 

    ### panel function for the inner nodes
    rval <- function(node) {
    
        node_vp <- viewport(x = unit(0.5, "npc"),
                        y = unit(0.5, "npc"),
                        width = unit(1, "strwidth", nstr) * 1.3, 
                        height = unit(3, "lines"),
 		        name = paste("node_inner", node$nodeID, sep = ""))
        pushViewport(node_vp)

        xell <- c(seq(0, 0.2, by = 0.01),
	          seq(0.2, 0.8, by = 0.05),
		  seq(0.8, 1, by = 0.01))
	yell <- sqrt(xell * (1-xell))
	
	lab <- getLabel1(node)
	fill <- rep(fill, length.out = 2)
	
        grid.polygon(x = unit(c(xell, rev(xell)), "npc"),
                     y = unit(c(yell, -yell)+0.5, "npc"),
                     gp = gpar(fill = fill[1]))
        grid.text(lab[1], y = unit(1.5 + 0.5 * pval, "lines"))
        if(pval) grid.text(lab[2], y = unit(1, "lines"))

        if (id) {
            nodeIDvp <- viewport(x = unit(0.5, "npc"), y = unit(1, "npc"),
	        width = max(unit(1, "lines"), unit(1.3, "strwidth", as.character(node$nodeID))),
	        height = max(unit(1, "lines"), unit(1.3, "strheight", as.character(node$nodeID))))
            pushViewport(nodeIDvp)
            grid.rect(gp = gpar(fill = fill[2]))
            grid.text(node$nodeID)
            popViewport()
        }
        upViewport()
    }
    
    return(rval)
}
class(node_inner) <- "grapcon_generator"

node_surv <- function(ctreeobj,
	 	      ylines = 2,
		      id = TRUE, ...)
{
    survobj <- response(ctreeobj)[[1]]
    if (!("Surv" %in% class(survobj))) 
        stop(sQuote("ctreeobj"), " is not a survival tree")

    ### panel function for Kaplan-Meier curves in nodes
    rval <- function(node) {
        km <- mysurvfit(survobj, weights = node$weights, ...)

        a <- dostep(km$time, km$surv)

        yscale <- c(0,1)
        xscale <- c(0, max(survobj[,1]))

        top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                           widths = unit(c(ylines, 1, 1), 
                                         c("lines", "null", "lines")),  
                           heights = unit(c(1, 1), c("lines", "null"))),
                           width = unit(1, "npc"), 
                           height = unit(1, "npc") - unit(2, "lines"),
			   name = paste("node_surv", node$nodeID, sep = ""))

        pushViewport(top_vp)
        grid.rect(gp = gpar(fill = "white", col = 0))

        ## main title
        top <- viewport(layout.pos.col=2, layout.pos.row=1)
        pushViewport(top)
	mainlab <- paste(ifelse(id, paste("Node", node$nodeID, "(n = "), "n = "),
	                 sum(node$weights), ifelse(id, ")", ""), sep = "")
        grid.text(mainlab)
        popViewport()
	
        plot <- viewport(layout.pos.col=2, layout.pos.row=2,
                         xscale=xscale, yscale=yscale,
			 name = paste("node_surv", node$nodeID, "plot", 
                         sep = ""))

        pushViewport(plot)
        grid.lines(a$x/max(survobj[,1]), a$y)
        grid.xaxis()
        grid.yaxis()
        grid.rect(gp = gpar(fill = "transparent"))
        upViewport(2)
    }

    return(rval)
}
class(node_surv) <- "grapcon_generator"

node_barplot <- function(ctreeobj,
                         col = "black",
      		         fill = NULL,
			 beside = NULL,
		         ymax = NULL,
		         ylines = NULL,
		         widths = 1,
		         gap = NULL,
			 reverse = NULL,
		         id = TRUE)
{   
    getMaxPred <- function(x) {
      mp <- max(x$prediction)
      mpl <- ifelse(x$terminal, 0, getMaxPred(x$left))
      mpr <- ifelse(x$terminal, 0, getMaxPred(x$right))
      return(max(c(mp, mpl, mpr)))
    }

    y <- response(ctreeobj)[[1]]
    
    if(is.factor(y) || class(y) == "was_ordered") {
        ylevels <- levels(y)
	if(is.null(beside)) beside <- if(length(ylevels) < 3) FALSE else TRUE
        if(is.null(ymax)) ymax <- if(beside) 1.1 else 1
	if(is.null(gap)) gap <- if(beside) 0.1 else 0
    } else {
        if(is.null(beside)) beside <- FALSE
        if(is.null(ymax)) ymax <- getMaxPred(ctreeobj@tree) * 1.1
        ylevels <- seq(along = ctreeobj@tree$prediction)
        if(length(ylevels) < 2) ylevels <- ""
	if(is.null(gap)) gap <- 1
    }
    if(is.null(reverse)) reverse <- !beside
    if(is.null(fill)) fill <- gray.colors(length(ylevels))
    if(is.null(ylines)) ylines <- if(beside) c(3, 2) else c(1.5, 2.5)

    ### panel function for barplots in nodes
    rval <- function(node) {
    
        ## parameter setup
        pred <- node$prediction
	if(reverse) {
	  pred <- rev(pred)
	  ylevels <- rev(ylevels)
	}
        np <- length(pred)
	nc <- if(beside) np else 1

	fill <- rep(fill, length.out = np)	
        widths <- rep(widths, length.out = nc)
	col <- rep(col, length.out = nc)
	ylines <- rep(ylines, length.out = 2)

	gap <- gap * sum(widths)
        yscale <- c(0, ymax)
        xscale <- c(0, sum(widths) + (nc+1)*gap)

        top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                           widths = unit(c(ylines[1], 1, ylines[2]), c("lines", "null", "lines")),
                           heights = unit(c(1, 1), c("lines", "null"))),
                           width = unit(1, "npc"), 
                           height = unit(1, "npc") - unit(2, "lines"),
			   name = paste("node_barplot", node$nodeID, sep = ""))

        pushViewport(top_vp)
        grid.rect(gp = gpar(fill = "white", col = 0))

        ## main title
        top <- viewport(layout.pos.col=2, layout.pos.row=1)
        pushViewport(top)
	mainlab <- paste(ifelse(id, paste("Node", node$nodeID, "(n = "), "n = "),
	                 sum(node$weights), ifelse(id, ")", ""), sep = "")
        grid.text(mainlab)
        popViewport()
	
        plot <- viewport(layout.pos.col=2, layout.pos.row=2,
                         xscale=xscale, yscale=yscale,
			 name = paste("node_barplot", node$nodeID, "plot", 
                         sep = ""))

        pushViewport(plot)
	
	if(beside) {
  	  xcenter <- cumsum(widths+gap) - widths/2
	  for (i in 1:np) {
            grid.rect(x = xcenter[i], y = 0, height = pred[i], 
                      width = widths[i],
	              just = c("center", "bottom"), default.units = "native",
	              gp = gpar(col = col[i], fill = fill[i]))
	  }
          if(length(xcenter) > 1) grid.xaxis(at = xcenter, label = FALSE)
	  grid.text(ylevels, x = xcenter, y = unit(-1, "lines"), 
                    just = c("center", "top"),
	            default.units = "native", check.overlap = TRUE)
          grid.yaxis()
	} else {
  	  ycenter <- cumsum(pred) - pred

	  for (i in 1:np) {
            grid.rect(x = xscale[2]/2, y = ycenter[i], height = min(pred[i], ymax - ycenter[i]), 
                      width = widths[1],
	              just = c("center", "bottom"), default.units = "native",
	              gp = gpar(col = col[i], fill = fill[i]))
	  }
          if(np > 1) {
	    grid.text(ylevels[1], x = unit(-1, "lines"), y = 0,
                      just = c("left", "center"), rot = 90,
	              default.units = "native", check.overlap = TRUE)
	    grid.text(ylevels[np], x = unit(-1, "lines"), y = ymax,
                      just = c("right", "center"), rot = 90,
	              default.units = "native", check.overlap = TRUE)
	  }
          if(np > 2) {
	    grid.text(ylevels[-c(1,np)], x = unit(-1, "lines"), y = ycenter[-c(1,np)],
                      just = "center", rot = 90,
	              default.units = "native", check.overlap = TRUE)
	  }
          grid.yaxis(main = FALSE)	
	}
	
        grid.rect(gp = gpar(fill = "transparent"))
        upViewport(2)
    }
    
    return(rval)
}
class(node_barplot) <- "grapcon_generator"

node_boxplot <- function(ctreeobj,
                         col = "black",
		         fill = "lightgray",
		         width = 0.5,
		         yscale = NULL,
		         ylines = 3,
			 cex = 0.5,
		         id = TRUE)
{
    y <- response(ctreeobj)[[1]]
    if (!is.numeric(y))
        stop(sQuote("ctreeobj"), " is not a regression tree")
    if (is.null(yscale)) 
        yscale <- range(y) + c(-0.1, 0.1) * diff(range(y))
         
    ### panel function for boxplots in nodes
    rval <- function(node) {
    
        ## parameter setup
	x <- boxplot(rep.int(y, node$weights), plot = FALSE)

        top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                           widths = unit(c(ylines, 1, 1), 
                                         c("lines", "null", "lines")),  
                           heights = unit(c(1, 1), c("lines", "null"))),
                           width = unit(1, "npc"), 
                           height = unit(1, "npc") - unit(2, "lines"),
			   name = paste("node_boxplot", node$nodeID, sep = ""))

        pushViewport(top_vp)
        grid.rect(gp = gpar(fill = "white", col = 0))

        ## main title
        top <- viewport(layout.pos.col=2, layout.pos.row=1)
        pushViewport(top)
	mainlab <- paste(ifelse(id, paste("Node", node$nodeID, "(n = "), "n = "),
	                 sum(node$weights), ifelse(id, ")", ""), sep = "")
        grid.text(mainlab)
        popViewport()
	
        plot <- viewport(layout.pos.col = 2, layout.pos.row = 2,
                         xscale = c(0, 1), yscale = yscale,
			 name = paste("node_boxplot", node$nodeID, "plot", 
                         sep = ""))

        pushViewport(plot)
	
	xl <- 0.5 - width/4
	xr <- 0.5 + width/4

        ## box & whiskers
        grid.lines(unit(c(xl, xr), "npc"), 
                   unit(x$stats[1], "native"), gp = gpar(col = col))
        grid.lines(unit(0.5, "npc"), 
                   unit(x$stats[1:2], "native"), gp = gpar(col = col, lty = 2))
        grid.rect(unit(0.5, "npc"), unit(x$stats[2], "native"), 
                  width = unit(width, "npc"), height = unit(diff(x$stats[c(2, 4)]), "native"),
                  just = c("center", "bottom"), 
                  gp = gpar(col = col, fill = fill))
        grid.lines(unit(c(0.5 - width/2, 0.5+width/2), "npc"), 
                   unit(x$stats[3], "native"), gp = gpar(col = col, lwd = 2))
        grid.lines(unit(0.5, "npc"), unit(x$stats[4:5], "native"), 
                   gp = gpar(col = col, lty = 2))
        grid.lines(unit(c(xl, xr), "npc"), unit(x$stats[5], "native"), 
                   gp = gpar(col = col))

        ## outlier
        n <- length(x$out)
        if (n > 0) {
            index <- 1:n ## which(x$out > yscale[1] & x$out < yscale[2])
            if (length(index) > 0)
                grid.points(unit(rep.int(0.5, length(index)), "npc"), 
                            unit(x$out[index], "native"),
                            size = unit(cex, "char"), gp = gpar(col = col))
        }
	
        grid.yaxis()
        grid.rect(gp = gpar(fill = "transparent"))
        upViewport(2)
    }
    
    return(rval)
}
class(node_boxplot) <- "grapcon_generator"

node_hist <- function(ctreeobj,
                      col = "black",
		      fill = "lightgray",
		      freq = FALSE,
		      horizontal = TRUE,
		      xscale = NULL,
		      ymax = NULL,
		      ylines = 3,
		      id = TRUE,
		      ...)
{
    y <- response(ctreeobj)[[1]]
    if (!is.numeric(y))
        stop(sQuote("ctreeobj"), " is not a regression tree")
    y <- rep.int(y, ctreeobj@tree$weights)
    yhist <- hist(y, plot = FALSE, ...)
    if (is.null(xscale)) 
        xscale <- range(yhist$breaks) + 
                        c(-0.05, 0.05) * diff(range(yhist$breaks))
    if (is.null(ymax)) {
        if (is.null(ymax)) 
            ymax <- if (freq) 0.7 * max(yhist$counts) 
                else 2.5 * max(yhist$density)
    }
    yscale <- c(0, ymax)
    
    if (horizontal) {
        yyy <- xscale
        xscale <- yscale
        yscale <- yyy
    }
         
    ### panel function for histograms in nodes
    rval <- function(node) {
    
        ## parameter setup
	yhist <- hist(rep.int(y, node$weights), plot = FALSE, ...)

        top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                           widths = unit(c(ylines, 1, 1), 
                                         c("lines", "null", "lines")),  
                           heights = unit(c(1, 1), c("lines", "null"))),
                           width = unit(1, "npc"), 
                           height = unit(1, "npc") - unit(2, "lines"),
			   name = paste("node_hist", node$nodeID, sep = ""))

        pushViewport(top_vp)
        grid.rect(gp = gpar(fill = "white", col = 0))

        ## main title
        top <- viewport(layout.pos.col=2, layout.pos.row=1)
        pushViewport(top)
	mainlab <- paste(ifelse(id, paste("Node", node$nodeID, "(n = "), "n = "),
	                 sum(node$weights), ifelse(id, ")", ""), sep = "")
        grid.text(mainlab)
        popViewport()
	
        plot <- viewport(layout.pos.col = 2, layout.pos.row = 2,
                         xscale = xscale, yscale = yscale,
			 name = paste("node_hist", node$nodeID, "plot", 
                         sep = ""))

        pushViewport(plot)
	
        ## histogram
        xpos <- yhist$breaks[-1]
	ypos <- 0
	yheight <- if (freq) yhist$counts else yhist$density
	xwidth <- diff(yhist$breaks)

        if (horizontal) {
              yyy <- xpos
              xpos <- ypos
              ypos <- yyy
              yyy <- xwidth
              xwidth <- -yheight
              yheight <- yyy
        }

	grid.rect(x = xpos, y = ypos,
	          width = xwidth, height = yheight,
		  just = c("right", "bottom"), default.units = "native",
		  gp = gpar(col = col, fill = fill))
	
        grid.xaxis()
        grid.yaxis()
        grid.rect(gp = gpar(fill = "transparent"))
        upViewport(2)
    }
    return(rval)
}
class(node_hist) <- "grapcon_generator"

node_density <- function(ctreeobj,
                         col = "black",
		         rug = TRUE,
		         horizontal = TRUE,
		         xscale = NULL,
		         yscale = NULL,
		         ylines = 3,
		         id = TRUE)
{
    y <- response(ctreeobj)[[1]]
    if (!is.numeric(y))
        stop(sQuote("ctreeobj"), " is not a regression tree")
    y <- rep.int(y, ctreeobj@tree$weights)
    ydens <- density(y)
    if (is.null(xscale)) 
        xscale <- range(ydens$x) + c(-0.05, 0.05) * diff(range(ydens$x))
    if (is.null(yscale)) {
        ymin <- if (rug) -max(ydens$y) * 0.1 else 0    
        yscale <- c(ymin, max(ydens$y) * 2.5)
    }

    xr <- xscale
    yr <- 0
    if (horizontal) {
        yyy <- xscale
        xscale <- yscale
        yscale <- yyy
    }
         
    ### panel function for density plots in nodes
    rval <- function(node) {
    
        ## parameter setup
	ydens <- density(rep.int(y, node$weights))

        top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                           widths = unit(c(ylines, 1, 1), 
                                         c("lines", "null", "lines")),  
                           heights = unit(c(1, 1), c("lines", "null"))),
                           width = unit(1, "npc"), 
                           height = unit(1, "npc") - unit(2, "lines"),
			   name = paste("node_density", node$nodeID, sep = ""))

        pushViewport(top_vp)
        grid.rect(gp = gpar(fill = "white", col = 0))

        ## main title
        top <- viewport(layout.pos.col=2, layout.pos.row=1)
        pushViewport(top)
	mainlab <- paste(ifelse(id, paste("Node", node$nodeID, "(n = "), "n = "),
	                 sum(node$weights), ifelse(id, ")", ""), sep = "")
        grid.text(mainlab)
        popViewport()
	
        plot <- viewport(layout.pos.col = 2, layout.pos.row = 2,
                         xscale = xscale, yscale = yscale,
			 name = paste("node_density", node$nodeID, "plot", 
                         sep = ""))

        pushViewport(plot)
	
        ## density
        yd <- ydens$y
	xd <- ydens$x

        if (horizontal) {
            yyy <- xd
            xd <- yd
            yd <- yyy
            yyy <- xr
            xr <- yr
            yr <- yyy
        }

	if (rug) {
            if (horizontal)
	        grid.rect(x = xscale[1], y = y[node$weights > 0],
	                  height = 0, width = xscale[1],
		          default.units = "native", 
                          just = c("right", "bottom"))
            else
                grid.rect(x = y[node$weights > 0], y = yscale[1],
	                  width = 0, height = abs(yscale[1]),
		          default.units = "native", 
                          just = c("center", "bottom"))

            grid.lines(x = xr, y = yr, gp = gpar(col = "lightgray"), 
                       default.units = "native")
            grid.lines(x = xr, y = yr, gp = gpar(col = "lightgray"), 
                       default.units = "native")
	}
	grid.lines(x = xd, y = yd, default.units = "native",
		   gp = gpar(col = col))
	
        grid.xaxis()
        grid.yaxis()
        grid.rect(gp = gpar(fill = "transparent"))
        upViewport(2)
    }
    
    return(rval)
}
class(node_density) <- "grapcon_generator"

node_terminal <- function(ctreeobj,
                          digits = 3,
		          abbreviate = FALSE,
		          fill = c("lightgray", "white"),
		          id = TRUE)
{
    getLabel1 <- function(x) {
        if (!x$terminal) return(rep.int("", 2))
        nlab <- paste("n =", sum(x$weights))
        ylab <- if (length(x$prediction) > 1)
                    paste("y =", paste("(", paste(round(x$prediction, digits),
	  	          collapse = ", "), ")", sep =""))
                else
	            paste("y =", round(x$prediction, digits))
      return(c(nlab, ylab))
    }

    maxstr <- function(node) {
        lab <- getLabel1(node)
        msl <- ifelse(node$terminal, "", maxstr(node$left))
        msr <- ifelse(node$terminal, "", maxstr(node$right))
        lab <- c(lab, msl, msr)
        return(lab[which.max(nchar(lab))])
    }

    nstr <- maxstr(ctreeobj@tree)

    ### panel function for simple n, Y terminal node labelling
    rval <- function(node) {
        fill <- rep(fill, length.out = 2)	

        node_vp <- viewport(x = unit(0.5, "npc"),   
                       y = unit(0.5, "npc"),   
                       width = unit(1, "strwidth", nstr) * 1.1,
                       height = unit(3, "lines"),
		       name = paste("node_terminal", node$nodeID, sep = ""))
        pushViewport(node_vp)

        lab <- getLabel1(node)
	
        grid.rect(gp = gpar(fill = fill[1]))
        grid.text(y = unit(2, "lines"), lab[1])
        grid.text(y = unit(1, "lines"), lab[2])

        if (id) {
            nodeIDvp <- viewport(x = unit(0.5, "npc"), y = unit(1, "npc"),
	        width = max(unit(1, "lines"), unit(1.3, "strwidth", as.character(node$nodeID))),
	        height = max(unit(1, "lines"), unit(1.3, "strheight", as.character(node$nodeID))))
            pushViewport(nodeIDvp)
            grid.rect(gp = gpar(fill = fill[2], lty = "solid"))
            grid.text(node$nodeID)
            popViewport()
	}
        upViewport()
    }
    return(rval)
}
class(node_terminal) <- "grapcon_generator"

edge_simple <- function(treeobj, digits = 3, abbreviate = FALSE)
{
    ### panel function for simple edge labelling
    function(split, ordered = FALSE, left = TRUE) {
  
        if (is.numeric(split)) 
            split <- round(split, digits = digits)
        if (is.character(split) & abbreviate > 0) 
            split <- abbreviate(split, as.numeric(abbreviate))

        if (!ordered) {
            if (length(split) > 1) 
                split <- paste("{", paste(split, collapse = ", "), 
    	                       "}", sep="")
        } else {
            ### <FIXME> phantom and . functions cannot be found by
            ###         codetools
            ### </FIXME>
            if (left) split <- as.expression(bquote(phantom(0) <= .(split)))
                else split <- as.expression(bquote(phantom(0) > .(split)))
        }
        grid.rect(gp = gpar(fill = "white", col = 0), 
                  width = unit(1, "strwidth", split)) 
        grid.text(split, just = "center")
    }
}
class(edge_simple) <- "grapcon_generator"

plotTree <- function(node, xlim, ylim, nx, ny, 
               terminal_panel, inner_panel, edge_panel,
	       tnex = 2, drop_terminal = TRUE, debug = FALSE) {

    ### the workhorse for plotting trees

    ### set up viewport for terminal node
    if (node$terminal) {
        x <- xlim[1] + diff(xlim)/2
        y <- ylim[1] + 0.5
       
        tn_vp <- viewport(x = unit(x, "native"),
                          y = unit(y, "native") - unit(0.5, "lines"),
                          width = unit(1, "native"), 
                          height = unit(tnex, "native") - unit(1, "lines"),
			  just = c("center", "top"),
                          name = paste("Node", node$nodeID, sep = ""))
        pushViewport(tn_vp)
        if (debug)
            grid.rect(gp = gpar(lty = "dotted", col = 4))
        terminal_panel(node) 
        upViewport()
        return(NULL)
    }    

    ### number of left leafs
    nl <- nterminal(node$left)

    ### number of right leafs
    nr <- nterminal(node$right)

    ### position of inner node
    x0 <- xlim[1] + (nl / (nl + nr)) * diff(xlim)
    y0 <- max(ylim)

    ### proportion of left terminal nodes in left node
    if (node$left$terminal) {
        lf <- 1/2
    } else {
        lf <- nterminal(node$left$left) / (nterminal(node$left$left) + 
                                           nterminal(node$left$right))
    }

    ### proportion of left terminal nodes in right node
    if (node$right$terminal) {
        rf <- 1/2
    } else {
        rf <- nterminal(node$right$left) / (nterminal(node$right$left) + 
                                            nterminal(node$right$right))
    }

    ### position of left and right daugher node
    x1l <- xlim[1] + (x0 - xlim[1]) * lf
    x1r <- x0 + (xlim[2] - x0) * rf
    
    if (!drop_terminal) {
        y1l <- y1r <- y0 - 1
    } else {
        y1l <- if (node$left$terminal) tnex - 0.5 else y0 - 1
        y1r <- if (node$right$terminal) tnex - 0.5 else y0 - 1
    }

    ### draw edges
    grid.lines(x = unit(c(x0, x1l), "native"), 
               y = unit(c(y0, y1l), "native"))
    grid.lines(x = unit(c(x0, x1r), "native"), 
               y = unit(c(y0, y1r), "native"))

    ### create viewport for inner node
    in_vp <- viewport(x = unit(x0, "native"),
                      y = unit(y0, "native"),
                      width = unit(1, "native"),
                      height = unit(1, "native") - unit(1, "lines"), 
                      name = paste("Node", node$nodeID, sep = ""))
    pushViewport(in_vp)
    if (debug)
        grid.rect(gp = gpar(lty = "dotted"))
    inner_panel(node)
    upViewport()

    ps <- node$psplit
    if (ps$ordered) {
        if (!is.null(attr(ps$splitpoint, "levels"))) {
            split <- attr(ps$splitpoint, "levels")[ps$splitpoint]
        } else {
            split <- ps$splitpoint
        }
    } else {
        ### <FIXME>: always to the left? </FIXME>
        split <- attr(ps$splitpoint, "levels")[as.logical(ps$splitpoint) & (ps$table > 0)]
    }


    ### position of labels
    y1lr <- max(y1l, y1r)
    ypos <- y0 - (y0 - y1lr) * 0.5
    xlpos <- x0 - (x0 - x1l) * 0.5 * (y0 - y1lr)/(y0 - y1l)
    xrpos <- x0 - (x0 - x1r) * 0.5 * (y0 - y1lr)/(y0 - y1r)

    ### setup left label
    lsp_vp <- viewport(x = unit(xlpos, "native"),
                       y = unit(ypos, "native"),
                       width = unit(xlpos - xrpos, "native"),
                       height = unit(1, "lines"), 
                       name =  paste("lEdge", node$nodeID, sep = ""))
    pushViewport(lsp_vp)
    if (debug)
        grid.rect(gp = gpar(lty = "dotted", col = 2))
    edge_panel(split, ordered = ps$ordered, left = TRUE)
    upViewport()

    ### setup right label
    if (ps$ordered) {
        if (!is.null(attr(ps$splitpoint, "levels"))) {
            split <- attr(ps$splitpoint, "levels")[ps$splitpoint]
        } else {
            split <- ps$splitpoint
        }
    } else {
        split <- attr(ps$splitpoint, "levels")[!as.logical(ps$splitpoint) & (ps$table > 0)]
    }

    rsp_vp <- viewport(x = unit(xrpos, "native"),
                       y = unit(ypos, "native"),
                       width = unit(xlpos - xrpos, "native"),
                       height = unit(1, "lines"),
                       name =  paste("rEdge", node$nodeID, sep = ""))
    pushViewport(rsp_vp) 
    if (debug)
        grid.rect(gp = gpar(lty = "dotted", col = 2))
    edge_panel(split, ordered = ps$ordered, left = FALSE)
    upViewport()

    plotTree(node$left, c(xlim[1], x0), c(y1l, 1), nx, ny, 
      terminal_panel, inner_panel, edge_panel,
      tnex = tnex, drop_terminal = drop_terminal, debug = debug)
    plotTree(node$right, c(x0, xlim[2]), c(y1r, 1), nx, ny,
      terminal_panel, inner_panel, edge_panel,
      tnex = tnex, drop_terminal = drop_terminal, debug = debug)
}


plot.BinaryTree <- function(x, main = NULL, type = c("extended", "simple"),
                            terminal_panel = NULL, tp_args = list(),
			    inner_panel = node_inner, ip_args = list(),
                            edge_panel = edge_simple, ep_args = list(),
			    drop_terminal = (type[1] == "extended"),
			    tnex = (type[1] == "extended") + 1, 
			    newpage = TRUE,
			    pop = TRUE,
			    ...) {

    ### plot BinaryTree objects

    ### extract tree
    ptr <- x@tree
    ### total number of terminal nodes
    nx <- nterminal(ptr)
    ### maximal depth of the tree
    ny <- maxdepth(ptr)

    ### compute default settings
    type <- match.arg(type)
    if (type == "simple") {
        if (is.null(terminal_panel)) 
            terminal_panel <- node_terminal
        if (is.null(tnex)) tnex <- 1
    } else {
        if (is.null(terminal_panel))
            terminal_panel <- switch(class(response(x)[[1]])[1],
	                             "Surv" = node_surv,
                                     "factor" = node_barplot,
                                     "was_ordered" = node_barplot,
                                     "ordered" = node_barplot,
                                     node_boxplot)
        if (is.null(tnex)) tnex <- 2
    }

    ## setup newpage
    if (newpage) grid.newpage()

    ## setup root viewport
    root_vp <- viewport(layout = grid.layout(3, 3, 
    			heights = unit(c(ifelse(is.null(main), 0, 3), 1, 1), 
                                      c("lines", "null", "lines")),
    			widths = unit(c(1, 1, 1), 
                                     c("lines", "null", "lines"))), 
    			name = "root")       
    pushViewport(root_vp)
  
    ## viewport for main title (if any)
    if (!is.null(main)) {
        main_vp <- viewport(layout.pos.col = 2, layout.pos.row = 1, 
                            name = "main")
        pushViewport(main_vp)
        grid.text(y = unit(1, "lines"), main, just = "center")
        upViewport()
    }

    ## setup viewport for tree
    tree_vp <- viewport(layout.pos.col = 2, layout.pos.row = 2, 
    			xscale = c(0, nx), yscale = c(0, ny + (tnex - 1)), 
                        name = "tree")
    pushViewport(tree_vp)

    ### setup panel functions (if necessary)
    ### the heuristic is as follows: If the first argument
    ### is `ctreeobj' than we assume a panel generating function, 
    ### otherwise the function is treated as a panel function
    if(inherits(terminal_panel, "grapcon_generator"))
      terminal_panel <- do.call("terminal_panel", c(list(x), as.list(tp_args)))
    if(inherits(inner_panel, "grapcon_generator"))
      inner_panel <- do.call("inner_panel", c(list(x), as.list(ip_args)))
    if(inherits(edge_panel, "grapcon_generator"))
      edge_panel <- do.call("edge_panel", c(list(x), as.list(ep_args)))


    if((nx <= 1 & ny <= 1)) {
      pushViewport(plotViewport(margins = rep(1.5, 4), name = paste("Node", ptr$nodeID, sep = "")))
      terminal_panel(ptr)    
    } else {
      ## call the workhorse
      plotTree(ptr,
        xlim = c(0, nx), ylim = c(0, ny - 0.5 + (tnex - 1)),
        nx = nx, ny = ny, 
        terminal_panel = terminal_panel,
        inner_panel = inner_panel,
        edge_panel = edge_panel,
        tnex = tnex,
        drop_terminal = drop_terminal,
        debug = FALSE)
    }
    upViewport()
    if (pop) popViewport() else upViewport()
}
