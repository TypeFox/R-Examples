.nobs_party <- function(party, id = 1L) {
  dat <- data_party(party, id = id)
  if("(weights)" %in% names(dat)) sum(dat[["(weights)"]]) else NROW(dat)
}

node_inner <- function(obj, id = TRUE, pval = TRUE, abbreviate = FALSE, fill = "white", gp = gpar())
{
  meta <- obj$data
  nam <- names(obj)

  extract_label <- function(node) {
    if(is.terminal(node)) return(rep.int("", 2L))

    varlab <- character_split(split_node(node), meta)$name
    if(abbreviate > 0L) varlab <- abbreviate(varlab, as.integer(abbreviate))

    ## FIXME: make more flexible rather than special-casing p-value
    if(pval) {
      pval <- suppressWarnings(try(!is.null(info_node(node)$p.value), silent = TRUE))
      pval <- if(inherits(pval, "try-error")) FALSE else pval
    }
    if(pval) {
      pvalue <- node$info$p.value
      plab <- ifelse(pvalue < 10^(-3L),
        paste("p <", 10^(-3L)),
    	paste("p =", round(pvalue, digits = 3L)))
    } else {
      plab <- ""
    }
    return(c(varlab, plab))
  }

  maxstr <- function(node) {
      lab <- extract_label(node)
      klab <- if(is.terminal(node)) "" else unlist(lapply(kids_node(node), maxstr))
      lab <- c(lab, klab)
      lab <- unlist(lapply(lab, function(x) strsplit(x, "\n")))
      lab <- lab[which.max(nchar(lab))]
      if(length(lab) < 1L) lab <- ""
      return(lab)
  }

  nstr <- maxstr(node_party(obj))
  if(nchar(nstr) < 6) nstr <- "aAAAAa"

  ### panel function for the inner nodes
  rval <- function(node) {  
    node_vp <- viewport(
      x = unit(0.5, "npc"),
      y = unit(0.5, "npc"),
      width = unit(1, "strwidth", nstr) * 1.3, 
      height = unit(3, "lines"),
      name = paste("node_inner", id_node(node), sep = ""),
      gp = gp
    )
    pushViewport(node_vp)

    xell <- c(seq(0, 0.2, by = 0.01),
  	      seq(0.2, 0.8, by = 0.05),
  	      seq(0.8, 1, by = 0.01))
    yell <- sqrt(xell * (1-xell))

    lab <- extract_label(node)
    fill <- rep(fill, length.out = 2L)

    grid.polygon(x = unit(c(xell, rev(xell)), "npc"),
        	 y = unit(c(yell, -yell)+0.5, "npc"),
        	 gp = gpar(fill = fill[1]))

    ## FIXME: something more general instead of pval ?
    grid.text(lab[1L], y = unit(1.5 + 0.5 * (lab[2L] != ""), "lines"))
    if(lab[2L] != "") grid.text(lab[2L], y = unit(1, "lines"))

    if(id) {
      nodeIDvp <- viewport(x = unit(0.5, "npc"), y = unit(1, "npc"),
        width = max(unit(1, "lines"), unit(1.3, "strwidth", nam[id_node(node)])),
        height = max(unit(1, "lines"), unit(1.3, "strheight", nam[id_node(node)])))
      pushViewport(nodeIDvp)
      grid.rect(gp = gpar(fill = fill[2]))
      grid.text(nam[id_node(node)])
      popViewport()
    }
    upViewport()
  }
  
  return(rval)
}
class(node_inner) <- "grapcon_generator"

node_terminal <- function(obj,
                          digits = 3,
		          abbreviate = FALSE,
		          fill = c("lightgray", "white"),
		          id = TRUE,
			  just = c("center", "top"),
			  top = 0.85,
			  align = c("center", "left", "right"),
			  gp = NULL,
			  FUN = NULL,
			  height = NULL,
			  width = NULL)
{
  nam <- names(obj)

  extract_label <- function(node) formatinfo_node(node, FUN = FUN, default = c("terminal", "node"))

  maxstr <- function(node) {
      lab <- extract_label(node)
      klab <- if(is.terminal(node)) "" else unlist(lapply(kids_node(node), maxstr))
      lab <- c(lab, klab)
      lab <- unlist(lapply(lab, function(x) strsplit(x, "\n")))
      return(lab[which.max(nchar(lab))])
  }

  nstr <- if(is.null(width)) maxstr(node_party(obj)) else paste(rep("a", width), collapse = "")

  just <- match.arg(just[1L], c("center", "centre", "top"))
  if(just == "centre") just <- "center"
  align <- match.arg(align[1L], c("center", "centre", "left", "right"))
  if(align == "centre") align <- "center"

  ### panel function for simple n, Y terminal node labeling
  rval <- function(node) {
    fill <- rep(fill, length.out = 2)	    

    lab <- extract_label(node)

    ## if gp is set, then an additional viewport may be
    ## required to appropriately evaluate strwidth unit
    if(!is.null(gp)) {
      outer_vp <- viewport(gp = gp)
      pushViewport(outer_vp)
    }
    
    if(is.null(height)) height <- length(lab) + 1L
    
    node_vp <- viewport(x = unit(0.5, "npc"),
      y = unit(if(just == "top") top else 0.5, "npc"),
      just = c("center", just),
      width = unit(1, "strwidth", nstr) * 1.1,
      height = unit(height, "lines"),
      name = paste("node_terminal", id_node(node), sep = ""),
      gp = if(is.null(gp)) gpar() else gp
    )
    pushViewport(node_vp)

    grid.rect(gp = gpar(fill = fill[1]))
      
    for(i in seq_along(lab)) grid.text(
      x = switch(align,
        "center" = unit(0.5, "npc"),
	"left"   = unit(1, "strwidth", "a"),
	"right"  = unit(1, "npc") - unit(1, "strwidth", "a")),
      y = unit(length(lab) - i + 1, "lines"), lab[i], just = align)

    if(id) {
      nodeIDvp <- viewport(x = unit(0.5, "npc"), y = unit(1, "npc"),
        width = max(unit(1, "lines"), unit(1.3, "strwidth", nam[id_node(node)])),
        height = max(unit(1, "lines"), unit(1.3, "strheight", nam[id_node(node)])))
      pushViewport(nodeIDvp)
      grid.rect(gp = gpar(fill = fill[2], lty = "solid"))
      grid.text(nam[id_node(node)])
      popViewport()
    }
    
    if(is.null(gp)) upViewport() else upViewport(2)
  }
  return(rval)
}
class(node_terminal) <- "grapcon_generator"

edge_simple <- function(obj, digits = 3, abbreviate = FALSE,
  justmin = Inf, just = c("alternate", "increasing", "decreasing", "equal"),
  fill = "white")
{
  meta <- obj$data

  justfun <- function(i, split) {
    myjust <- if(mean(nchar(split)) > justmin) {
      match.arg(just, c("alternate", "increasing", "decreasing", "equal"))
    } else {
      "equal"
    }
    k <- length(split)
    rval <- switch(myjust,
      "equal" = rep.int(0, k),
      "alternate" = rep(c(0.5, -0.5), length.out = k),
      "increasing" = seq(from = -k/2, to =  k/2, by = 1),
      "decreasing" = seq(from =  k/2, to = -k/2, by = -1)
    )
    unit(0.5, "npc") + unit(rval[i], "lines")
  }

  ### panel function for simple edge labelling
  function(node, i) {
    split <- character_split(split_node(node), meta, digits = digits)$levels
    y <- justfun(i, split)
    split <- split[i]
    # try() because the following won't work for split = "< 10 Euro", for example.
    if(any(grep(">", split) > 0) | any(grep("<", split) > 0)) {
      tr <- suppressWarnings(try(parse(text = paste("phantom(0)", split)), silent = TRUE))
      if(!inherits(tr, "try-error")) split <- tr
    }
    grid.rect(y = y, gp = gpar(fill = fill, col = 0), width = unit(1, "strwidth", split))
    grid.text(split, y = y, just = "center")
  }
}
class(edge_simple) <- "grapcon_generator"

.plot_node <- function(node, xlim, ylim, nx, ny, 
               terminal_panel, inner_panel, edge_panel,
	       tnex = 2, drop_terminal = TRUE, debug = FALSE) {

    ### the workhorse for plotting trees
 
    ### set up viewport for terminal node
    if (is.terminal(node)) {
        x <- xlim[1] + diff(xlim)/2
        y <- ylim[1] + 0.5
       
        tn_vp <- viewport(x = unit(x, "native"),
                          y = unit(y, "native") - unit(0.5, "lines"),
                          width = unit(1, "native"), 
                          height = unit(tnex, "native") - unit(1, "lines"),
			  just = c("center", "top"),
                          name = paste("Node", id_node(node), sep = ""))
        pushViewport(tn_vp)
        if (debug)
            grid.rect(gp = gpar(lty = "dotted", col = 4))
        terminal_panel(node) 
        upViewport()
        return(NULL)
    }    

    ## convenience function for computing relative position of splitting node
    pos_frac <- function(node) {
      if(is.terminal(node)) 0.5 else {
        width_kids <- sapply(kids_node(node), width)
        nk <- length(width_kids)
        rval <- if(nk %% 2 == 0) sum(width_kids[1:(nk/2)]) else
	  mean(cumsum(width_kids)[nk/2 + c(-0.5, 0.5)])
	rval/sum(width_kids)
      }
    }

    ## extract information
    split <- split_node(node)
    kids <- kids_node(node)
    width_kids <- sapply(kids, width)
    nk <- length(width_kids)

    ### position of inner node
    x0 <- xlim[1] + pos_frac(node) * diff(xlim)
    y0 <- max(ylim)

    ### relative positions of kids
    xfrac <- sapply(kids, pos_frac)
    x1lim <- xlim[1] + cumsum(c(0, width_kids))/sum(width_kids) * diff(xlim)
    x1 <- x1lim[1:nk] + xfrac * diff(x1lim)
    if (!drop_terminal) {
        y1 <- rep(y0 - 1, nk)
    } else {
        y1 <- ifelse(sapply(kids, is.terminal), tnex - 0.5, y0 - 1)
    }

    ### draw edges
    for(i in 1:nk) grid.lines(x = unit(c(x0, x1[i]), "native"), y = unit(c(y0, y1[i]), "native"))

    ### create viewport for inner node
    in_vp <- viewport(x = unit(x0, "native"),
                      y = unit(y0, "native"),
                      width = unit(1, "native"),
                      height = unit(1, "native") - unit(1, "lines"), 
                      name = paste("Node", id_node(node), sep = ""))
    pushViewport(in_vp)
    if(debug) grid.rect(gp = gpar(lty = "dotted"))
    inner_panel(node)
    upViewport()

    ### position of labels
    y1max <- max(y1)
    ypos <- y0 - (y0 - y1max) * 0.5
    xpos <- x0 - (x0 - x1) * 0.5 * (y0 - y1max)/(y0 - y1)

    ### setup labels
    for(i in 1:nk) {
      sp_vp <- viewport(x = unit(xpos[i], "native"),
                        y = unit(ypos, "native"),
                        width = unit(diff(x1lim)[i], "native"),
                        height = unit(1, "lines"), 
                        name =  paste("edge", id_node(node), "-", i, sep = ""))
      pushViewport(sp_vp)
      if(debug) grid.rect(gp = gpar(lty = "dotted", col = 2))
      edge_panel(node, i)
      upViewport()
    }

    ## call workhorse for kids
    for(i in 1:nk) .plot_node(kids[[i]],
      c(x1lim[i], x1lim[i+1]), c(y1[i], 1), nx, ny, 
      terminal_panel, inner_panel, edge_panel,
      tnex = tnex, drop_terminal = drop_terminal, debug = debug)
}


plot.party <- function(x, main = NULL,
                       terminal_panel = node_terminal, tp_args = list(),
		       inner_panel = node_inner, ip_args = list(),
                       edge_panel = edge_simple, ep_args = list(),
		       drop_terminal = FALSE, tnex = 1, 
		       newpage = TRUE, pop = TRUE, gp = gpar(), ...)
{

    ### extract tree
    node <- node_party(x)
    ### total number of terminal nodes
    nx <- width(node)
    ### maximal depth of the tree
    ny <- depth(node, root = TRUE)

    ## setup newpage
    if (newpage) grid.newpage()

    ## setup root viewport
    root_vp <- viewport(layout = grid.layout(3, 3, 
    			heights = unit(c(ifelse(is.null(main), 0, 3), 1, 1), 
                                      c("lines", "null", "lines")),
    			widths = unit(c(1, 1, 1), 
                                     c("lines", "null", "lines"))), 
    			name = "root",
			gp = gp)       
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
    if(inherits(terminal_panel, "grapcon_generator"))
      terminal_panel <- do.call("terminal_panel", c(list(x), as.list(tp_args)))
    if(inherits(inner_panel, "grapcon_generator"))
      inner_panel <- do.call("inner_panel", c(list(x), as.list(ip_args)))
    if(inherits(edge_panel, "grapcon_generator"))
      edge_panel <- do.call("edge_panel", c(list(x), as.list(ep_args)))


    if((nx <= 1 & ny <= 1)) {
      pushViewport(plotViewport(margins = rep(1.5, 4), name = paste("Node", id_node(node), sep = "")))
      terminal_panel(node)
    } else {
      ## call the workhorse
      .plot_node(node,
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

plot.constparty <- function(x, main = NULL,
                        terminal_panel = NULL, tp_args = list(),
	 	        inner_panel = node_inner, ip_args = list(),
                        edge_panel = edge_simple, ep_args = list(),
		        type = c("extended", "simple"), drop_terminal = NULL, tnex = NULL, 
		        newpage = TRUE, pop = TRUE, gp = gpar(), ...)
{
    ### compute default settings
    type <- match.arg(type)
    if (type == "simple") {
	x <- as.simpleparty(x)
        if (is.null(terminal_panel)) 
            terminal_panel <- node_terminal
        if (is.null(tnex)) tnex <- 1
        if (is.null(drop_terminal)) drop_terminal <- FALSE
        if (is.null(tp_args) || length(tp_args) < 1L) tp_args <- list(
	  FUN = .make_formatinfo_simpleparty(x, digits = getOption("digits") - 4L, sep = "\n"))
    } else {
        if (is.null(terminal_panel)) {
	    cl <- class(x$fitted[["(response)"]])
	    if("factor" %in% cl) {
	        terminal_panel <- node_barplot 
	    } else if("Surv" %in% cl) {
	        terminal_panel <- node_surv
            } else if ("data.frame" %in% cl) {
                terminal_panel <- node_mvar
                if (is.null(tnex)) tnex <- 2 * NCOL(x$fitted[["(response)"]])
            } else {
	        terminal_panel <- node_boxplot
	    }
	}
        if (is.null(tnex)) tnex <- 2
        if (is.null(drop_terminal)) drop_terminal <- TRUE
    }

    plot.party(x, main = main,
      terminal_panel = terminal_panel, tp_args = tp_args,
      inner_panel = inner_panel, ip_args = ip_args,
      edge_panel = edge_panel, ep_args = ep_args,
      drop_terminal = drop_terminal, tnex = tnex,
      newpage = newpage, pop = pop, gp = gp, ...)
}

node_barplot <- function(obj,
                         col = "black",
      		         fill = NULL,
			 bg = "white",
			 beside = NULL,
		         ymax = NULL,
		         ylines = NULL,
		         widths = 1,
		         gap = NULL,
			 reverse = NULL,
		         id = TRUE,
                         mainlab = NULL,
			 gp = gpar())
{   
    ## extract response
    y <- obj$fitted[["(response)"]]
    stopifnot(is.factor(y) || isTRUE(all.equal(round(y), y)) || is.data.frame(y))
    
    ## FIXME: This could be avoided by
    ##   predict_party(obj, nodeids(obj, terminal = TRUE), type = "prob")
    ## but only for terminal nodes                  ^^^^
    probs_and_n <- function(x) {
      y1 <- x$fitted[["(response)"]]
      if(!is.factor(y1)) {
        if(is.data.frame(y1)) {
	  y1 <- t(as.matrix(y1))
	} else {
          y1 <- factor(y1, levels = min(y):max(y))
	}
      }
      w <- x$fitted[["(weights)"]]
      if(is.null(w)) w <- rep.int(1L, length(y1))
      sumw <- if(is.factor(y1)) tapply(w, y1, sum) else drop(y1 %*% w)
      sumw[is.na(sumw)] <- 0
      prob <- c(sumw/sum(w), sum(w))
      names(prob) <- c(if(is.factor(y1)) levels(y1) else rownames(y1), "nobs")
      prob
    }
    probs <- do.call("rbind", nodeapply(obj, nodeids(obj), probs_and_n, by_node = FALSE))
    nobs <- probs[, "nobs"]
    probs <- probs[, -ncol(probs), drop = FALSE]
    
    if(is.factor(y)) {
        ylevels <- levels(y)
	if(is.null(beside)) beside <- if(length(ylevels) < 3L) FALSE else TRUE
        if(is.null(ymax)) ymax <- if(beside) 1.1 else 1
	if(is.null(gap)) gap <- if(beside) 0.1 else 0
    } else {
        if(is.null(beside)) beside <- TRUE
        if(is.null(ymax)) ymax <- if(beside) max(probs) * 1.1 else max(probs)
        ylevels <- colnames(probs)
        if(length(ylevels) < 2) ylevels <- ""
	if(is.null(gap)) gap <- if(beside) 0.1 else 0
    }
    if(is.null(reverse)) reverse <- !beside
    if(is.null(fill)) fill <- gray.colors(length(ylevels))
    if(is.null(ylines)) ylines <- if(beside) c(3, 2) else c(1.5, 2.5)

    ### panel function for barplots in nodes
    rval <- function(node) {
    
        ## id
	nid <- id_node(node)
    
        ## parameter setup
        pred <- probs[nid,]
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
			   name = paste0("node_barplot", nid),
			   gp = gp)

        pushViewport(top_vp)
        grid.rect(gp = gpar(fill = bg, col = 0))

        ## main title
        top <- viewport(layout.pos.col=2, layout.pos.row=1)
        pushViewport(top)
        if (is.null(mainlab)) {	
	  mainlab <- if(id) {
	    function(id, nobs) sprintf("Node %s (n = %s)", id, nobs)
  	  } else {
	    function(id, nobs) sprintf("n = %s", nobs)
	  }
        }
	if (is.function(mainlab)) {
          mainlab <- mainlab(names(obj)[nid], nobs[nid])
	}
        grid.text(mainlab)
        popViewport()
	
        plot <- viewport(layout.pos.col=2, layout.pos.row=2,
                         xscale=xscale, yscale=yscale,
			 name = paste0("node_barplot", node$nodeID, "plot"),
			 clip = FALSE)

        pushViewport(plot)
	
	if(beside) {
  	  xcenter <- cumsum(widths+gap) - widths/2
          if(length(xcenter) > 1) grid.xaxis(at = xcenter, label = FALSE)
	  grid.text(ylevels, x = xcenter, y = unit(-1, "lines"), 
                    just = c("center", "top"),
	            default.units = "native", check.overlap = TRUE)
          grid.yaxis()
          grid.rect(gp = gpar(fill = "transparent"))
	  grid.clip()
	  for (i in 1:np) {
            grid.rect(x = xcenter[i], y = 0, height = pred[i], 
                      width = widths[i],
	              just = c("center", "bottom"), default.units = "native",
	              gp = gpar(col = col[i], fill = fill[i]))
	  }
	} else {
  	  ycenter <- cumsum(pred) - pred

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

          grid.clip()
          grid.rect(gp = gpar(fill = "transparent"))
	  for (i in 1:np) {
            grid.rect(x = xscale[2]/2, y = ycenter[i], height = min(pred[i], ymax - ycenter[i]), 
                      width = widths[1],
	              just = c("center", "bottom"), default.units = "native",
	              gp = gpar(col = col[i], fill = fill[i]))
	  }
	}
	grid.rect(gp = gpar(fill = "transparent"))

	
        upViewport(2)
    }
    
    return(rval)
}
class(node_barplot) <- "grapcon_generator"

node_boxplot <- function(obj,
                         col = "black",
		         fill = "lightgray",
			 bg = "white",
		         width = 0.5,
		         yscale = NULL,
		         ylines = 3,
			 cex = 0.5,
		         id = TRUE,
                         mainlab = NULL, 
			 gp = gpar())
{
    y <- obj$fitted[["(response)"]]
    stopifnot(is.numeric(y))

    if (is.null(yscale)) 
        yscale <- range(y) + c(-0.1, 0.1) * diff(range(y))
         
    ### panel function for boxplots in nodes
    rval <- function(node) {

        ## extract data
	nid <- id_node(node)
	dat <- data_party(obj, nid)
	yn <- dat[["(response)"]]
	wn <- dat[["(weights)"]]
	if(is.null(wn)) wn <- rep(1, length(yn))
    
        ## parameter setup
	x <- boxplot(rep.int(yn, wn), plot = FALSE)

        top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                           widths = unit(c(ylines, 1, 1), 
                                         c("lines", "null", "lines")),  
                           heights = unit(c(1, 1), c("lines", "null"))),
                           width = unit(1, "npc"), 
                           height = unit(1, "npc") - unit(2, "lines"),
			   name = paste("node_boxplot", nid, sep = ""),
			   gp = gp)

        pushViewport(top_vp)
        grid.rect(gp = gpar(fill = bg, col = 0))

        ## main title
        top <- viewport(layout.pos.col=2, layout.pos.row=1)
        pushViewport(top)
        if (is.null(mainlab)) {	
	  mainlab <- if(id) {
	    function(id, nobs) sprintf("Node %s (n = %s)", id, nobs)
  	  } else {
	    function(id, nobs) sprintf("n = %s", nobs)
	  }
        }
	if (is.function(mainlab)) {
          mainlab <- mainlab(names(obj)[nid], sum(wn))
	}
        grid.text(mainlab)
        popViewport()
	
        plot <- viewport(layout.pos.col = 2, layout.pos.row = 2,
                         xscale = c(0, 1), yscale = yscale,
			 name = paste0("node_boxplot", nid, "plot"),
			 clip = FALSE)

        pushViewport(plot)
	
        grid.yaxis()
        grid.rect(gp = gpar(fill = "transparent"))
	grid.clip()

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
	
        upViewport(2)
    }
    
    return(rval)
}
class(node_boxplot) <- "grapcon_generator"

node_surv <- function(obj, col = "black", bg = "white", ylines = 2,
		      id = TRUE, mainlab = NULL, gp = gpar(), ...)
{
    ## extract response
    y <- obj$fitted[["(response)"]]
    stopifnot(inherits(y, "Surv"))

    ## helper functions
    mysurvfit <- function(y, weights, ...)
        survfit(y ~ 1, weights = weights) 
        ### structure(
        ###   survival:::survfitKM(x = gl(1, NROW(y)), y = y, casewt = weights, ...),
	### class = "survfit")

    dostep <- function(x, y) {
        ### create a step function based on x, y coordinates
        ### modified from `survival:print.survfit'
        if (is.na(x[1] + y[1])) {
            x <- x[-1]
            y <- y[-1]
        }
        n <- length(x)
        if (n > 2) {  
            # replace verbose horizonal sequences like
            # (1, .2), (1.4, .2), (1.8, .2), (2.3, .2), (2.9, .2), (3, .1)
            # with (1, .2), (3, .1).  They are slow, and can smear the looks
            # of the line type.
            dupy <- c(TRUE, diff(y[-n]) !=0, TRUE)
            n2 <- sum(dupy)

            #create a step function
            xrep <- rep(x[dupy], c(1, rep(2, n2-1)))
            yrep <- rep(y[dupy], c(rep(2, n2-1), 1))
            RET <- list(x = xrep, y = yrep)
        } else {
            if (n == 1) {
                RET <- list(x = x, y = y)
            } else {
                RET <- list(x = x[c(1,2,2)], y = y[c(1,1,2)])
            }
        }
        return(RET)
    }

    ### panel function for Kaplan-Meier curves in nodes
    rval <- function(node) {

        ## extract data
	nid <- id_node(node)
	dat <- data_party(obj, nid)
	yn <- dat[["(response)"]]
	wn <- dat[["(weights)"]]
	if(is.null(wn)) wn <- rep(1, NROW(yn))

        ## get Kaplan-Meier curver in node
        km <- mysurvfit(yn, weights = wn, ...)
        a <- dostep(km$time, km$surv)

        ## set up plot
        yscale <- c(0, 1)
        xscale <- c(0, max(y[,1]))

        top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                           widths = unit(c(ylines, 1, 1), 
                                         c("lines", "null", "lines")),  
                           heights = unit(c(1, 1), c("lines", "null"))),
                           width = unit(1, "npc"), 
                           height = unit(1, "npc") - unit(2, "lines"),
			   name = paste("node_surv", nid, sep = ""), gp = gp)

        pushViewport(top_vp)
        grid.rect(gp = gpar(fill = bg, col = 0))

        ## main title
        top <- viewport(layout.pos.col=2, layout.pos.row=1)
        pushViewport(top)
        if (is.null(mainlab)) {	
	  mainlab <- if(id) {
	    function(id, nobs) sprintf("Node %s (n = %s)", id, nobs)
  	  } else {
	    function(id, nobs) sprintf("n = %s", nobs)
	  }
        }
	if (is.function(mainlab)) {
          mainlab <- mainlab(nid, sum(wn))
	}
        grid.text(mainlab)
        popViewport()
	
        plot <- viewport(layout.pos.col=2, layout.pos.row=2,
                         xscale=xscale, yscale=yscale,
			 name = paste0("node_surv", nid, "plot"),
			 clip = FALSE)

        pushViewport(plot)
        grid.xaxis()
        grid.yaxis()
        grid.rect(gp = gpar(fill = "transparent"))
	grid.clip()
        grid.lines(a$x/max(a$x), a$y, gp = gpar(col = col))
        upViewport(2)
    }

    return(rval)
}
class(node_surv) <- "grapcon_generator"

node_ecdf <- function(obj, col = "black", bg = "white", ylines = 2,
		      id = TRUE, mainlab = NULL, gp = gpar(), ...)
{
    ## extract response
    y <- obj$fitted[["(response)"]]
    stopifnot(inherits(y, "numeric") || inherits(y, "integer"))

    dostep <- function(f) {
        x <- knots(f)
        y <- f(x)
        ### create a step function based on x, y coordinates
        ### modified from `survival:print.survfit'
        if (is.na(x[1] + y[1])) {
            x <- x[-1]
            y <- y[-1]
        }
        n <- length(x)
        if (n > 2) {  
            # replace verbose horizonal sequences like
            # (1, .2), (1.4, .2), (1.8, .2), (2.3, .2), (2.9, .2), (3, .1)
            # with (1, .2), (3, .1).  They are slow, and can smear the looks
            # of the line type.
            dupy <- c(TRUE, diff(y[-n]) !=0, TRUE)
            n2 <- sum(dupy)

            #create a step function
            xrep <- rep(x[dupy], c(1, rep(2, n2-1)))
            yrep <- rep(y[dupy], c(rep(2, n2-1), 1))
            RET <- list(x = xrep, y = yrep)
        } else {
            if (n == 1) {
                RET <- list(x = x, y = y)
            } else {
                RET <- list(x = x[c(1,2,2)], y = y[c(1,1,2)])
            }
        }
        return(RET)
    }

    ### panel function for ecdf in nodes
    rval <- function(node) {

        ## extract data
	nid <- id_node(node)
	dat <- data_party(obj, nid)
	yn <- dat[["(response)"]]
	wn <- dat[["(weights)"]]
	if(is.null(wn)) wn <- rep(1, NROW(yn))

        ## get ecdf in node
        f <- .pred_ecdf(yn, wn)
        a <- dostep(f)

        ## set up plot
        yscale <- c(0, 1)
        xscale <- range(y)
        a$x <- c(xscale[1], a$x[1], a$x, xscale[2])
        a$x <- a$x - min(a$x)
        a$x <- a$x / max(a$x)
        a$y <- c(0, 0, a$y, 1)

        top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                           widths = unit(c(ylines, 1, 1), 
                                         c("lines", "null", "lines")),  
                           heights = unit(c(1, 1), c("lines", "null"))),
                           width = unit(1, "npc"), 
                           height = unit(1, "npc") - unit(2, "lines"),
			   name = paste("node_ecdf", nid, sep = ""), gp = gp)

        pushViewport(top_vp)
        grid.rect(gp = gpar(fill = bg, col = 0))

        ## main title
        top <- viewport(layout.pos.col=2, layout.pos.row=1)
        pushViewport(top)
        if (is.null(mainlab)) {	
	  mainlab <- if(id) {
	    function(id, nobs) sprintf("Node %s (n = %s)", id, nobs)
  	  } else {
	    function(id, nobs) sprintf("n = %s", nobs)
	  }
        }
	if (is.function(mainlab)) {
          mainlab <- mainlab(nid, sum(wn))
	}
        grid.text(mainlab)
        popViewport()
	
        plot <- viewport(layout.pos.col=2, layout.pos.row=2,
                         xscale=xscale, yscale=yscale,
			 name = paste0("node_surv", nid, "plot"),
			 clip = FALSE)

        pushViewport(plot)
        grid.xaxis()
        grid.yaxis()
        grid.rect(gp = gpar(fill = "transparent"))
	grid.clip()
        grid.lines(a$x, a$y, gp = gpar(col = col))
        upViewport(2)
    }

    return(rval)
}
class(node_ecdf) <- "grapcon_generator"



node_mvar <- function(obj, which = NULL, id = TRUE, pop = TRUE, ylines = NULL, mainlab = NULL, varlab = TRUE, bg = "white", ...)
{  
  ## obtain dependent variables
  y <- obj$fitted[["(response)"]]

  ## fitted node ids
  fitted <- obj$fitted[["(fitted)"]]

  ## number of panels needed
  if(is.null(which)) which <- 1L:NCOL(y)
  k <- length(which)

  rval <- function(node) {
    
    tid <- id_node(node)
    nobs <- .nobs_party(obj, id = tid)

    ## set up top viewport
    top_vp <- viewport(layout = grid.layout(nrow = k, ncol = 2,
		       widths = unit(c(ylines, 1), c("lines", "null")), heights = unit(k, "null")),
		       width = unit(1, "npc"), height = unit(1, "npc") - unit(2, "lines"),
		       name = paste("node_mvar", tid, sep = ""))
    pushViewport(top_vp)
    grid.rect(gp = gpar(fill = bg, col = 0))

    ## main title
    if (is.null(mainlab)) { 
      mainlab <- if(id) {
    	function(id, nobs) sprintf("Node %s (n = %s)", id, nobs)
      } else {
    	function(id, nobs) sprintf("n = %s", nobs)
      }
    }
    if (is.function(mainlab)) {
      mainlab <- mainlab(tid, nobs)
    }
    
    for(i in 1L:k) {
      tmp <- obj
      tmp$fitted[["(response)"]] <- y[,which[i]]
      if(varlab) {
        nm <- names(y)[which[i]]
        if(i == 1L) nm <- paste(mainlab, nm, sep = ": ")
      } else {
        nm <- if(i == 1L) mainlab else ""
      }
      pfun <- switch(sapply(y, class)[which[i]],
                     "Surv" = node_surv(tmp, id = id, mainlab = nm, ...),
                     "factor" = node_barplot(tmp, id = id, mainlab = nm,  ...),
                     "ordered" = node_barplot(tmp, id = id, mainlab = nm, ...),
                     node_boxplot(tmp, id = id, mainlab = nm, ...))
      ## select panel
      plot_vpi <- viewport(layout.pos.col = 2L, layout.pos.row = i)
      pushViewport(plot_vpi)

      ## call panel function
      pfun(node)

      if(pop) popViewport() else upViewport()
    }
    if(pop) popViewport() else upViewport()
  }
  
  return(rval)
}
class(node_mvar) <- "grapcon_generator"
