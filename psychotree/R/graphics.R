##### panel generating graphic functions for various tree models

## profile plot visualization function
node_profileplot <- function(mobobj, what = c("items", "thresholds", "discriminations"),
  parg = list(type = NULL, ref = NULL, alias = TRUE), id = TRUE, names = FALSE,
  abbreviate = TRUE, index = TRUE, ref = TRUE, col = "black", border = col,
  linecol = "black", refcol = "lightgray", cex = 0.5, pch = 21, xscale = NULL, yscale = NULL, ylines = 2, ...)
{
  ## check input
  what <- match.arg(what)
  if (what == "thresholds") type <- parg$type
  refpar <- parg$ref
  alias <- if (is.null(parg$alias)) TRUE else parg$alias
  addargs <- list(...)
  if ("worth" %in% names(addargs)) warning("the argument 'worth' is deprecated and not longer used")

  ## node ids
  node <- nodeids(mobobj, terminal = FALSE)
  
  ## get all coefficients 
  if (what == "items") {
    cf <- apply_to_models(mobobj, node, FUN = function(z) coef(itempar(z, ref = refpar, alias = alias, vcov = FALSE)))
  } else if (what == "thresholds") {
    cf <- apply_to_models(mobobj, node, FUN = function(z) coef(threshpar(z, type = type, ref = refpar, alias = alias, vcov = FALSE), type = "matrix"))
  } else {
    cf <- apply_to_models(mobobj, node, FUN = function(z) coef(discrpar(z, ref = refpar, alias = alias, vcov = FALSE)))
  }
  names(cf) <- node

  ## labeling
  if (isTRUE(names)) {
    nms <- if (what != "thresholds") lapply(cf, names) else lapply(cf, rownames)
  } else if (is.character(names)) {
    nms <- split(rep(names, length(node)), f = rep(1:length(node), each = length(names)))
  } else {
    ncf <- lapply(cf, NROW)
    nms <- lapply(ncf, function(m) {
      lab <- rep("", m)
      lab[c(1, m)] <- c(1, m)
      pr <- pretty(1:m, n = 4)
      pr <- pr[pr > 1 & pr < m]
      lab[pr] <- pr    
      lab
    })
    abbreviate <- FALSE
  }
  
  ## abbreviation
  if (is.logical(abbreviate)) {
    nlab <- max(unlist(lapply(nms, function (j) nchar(j))))
    abbreviate <- if (abbreviate) as.numeric(cut(nlab, c(-Inf, 1.5, 4.5, 7.5, Inf))) else nlab
  }
  nms <- lapply(nms, function (j) abbreviate(j, abbreviate))

  ## axis scale
  if (index) {
    x <- if (what == "thresholds") 1:nrow(cf[[1]]) else 1:length(cf[[1]])
    if (is.null(xscale)) xscale <- range(x) + c(-0.1, 0.1) * diff(range(x))
  } else {
    if (what == "thresholds") {
      x <- 0:(ncol(cf[[1]]) - 1)
      if (is.null(xscale)) xscale <- c(-1, ncol(cf[[1]]))
    } else {
      x <- rep(0, length(cf[[1]]))
      if (is.null(xscale)) xscale <- c(-1, 1)
    }
  }
  rg <- range(unlist(cf)[is.finite(unlist(cf))], na.rm = TRUE)
  r <- diff(rg)
  if (!r) r <- 1
  if (is.null(yscale)) yscale <- rg + c(-0.1, 0.1) * r

  ## panel function for profile plots in nodes
  panelfun <- function (node) {

    ## node index
    idn <- id_node(node)
    
    ## get cfs and labels
    cfi <- cf[[idn]]
    if(any(!is.finite(cfi))) {
      cfi[cfi < 0 & !is.finite(cfi)] <- yscale[1]
      cfi[cfi > 0 & !is.finite(cfi)] <- yscale[2]
    }
    nmsi <- nms[[idn]]

    ## viewport setup
    top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                           widths = unit(c(ylines, 1, 1), c("lines", "null", "lines")),  
                           heights = unit(c(1, 1), c("lines", "null"))),
                       width = unit(1, "npc"), height = unit(1, "npc") - unit(2, "lines"),
                       name = paste("node_profileplot", idn, sep = ""))
    pushViewport(top_vp)
    grid.rect(gp = gpar(fill = "white", col = 0))

    ## main title
    top <- viewport(layout.pos.col = 2, layout.pos.row = 1)
    pushViewport(top)
    mainlab <- paste(ifelse(id, paste("Node", idn, "(n = "), ""),
                     info_node(node)$nobs, ifelse(idn, ")", ""), sep = "")
    grid.text(mainlab)
    popViewport()

    ## actual plot  
    plot_vpi <- viewport(layout.pos.col = 2, layout.pos.row = 2, xscale = xscale, yscale = yscale, 
                               name = paste("node_profileplot", idn, "plot", sep = ""))
    pushViewport(plot_vpi)
    grid.lines(xscale, c(mean(cfi), mean(cfi)), gp = gpar(col = refcol), default.units = "native")
    if(index) {
      if (what == "thresholds") {
        for (j in 1:ncol(cfi)) {
          grid.lines(x, cfi[, j], gp = gpar(col = linecol, lty = 2), default.units = "native")
          grid.text(label = paste0("C", j), x, cfi[, j], gp = gpar(col = col), default.units = "native")
        }
      } else {
        grid.lines(x, cfi, gp = gpar(col = linecol, lty = 2), default.units = "native")
        grid.points(x, cfi, gp = gpar(col = border, fill = col, cex = cex), pch = pch, default.units = "native")
      }
      grid.xaxis(at = x, label = nmsi)
    } else {	
      if (what == "thresholds") {
        for (j in 1:ncol(cfi)) grid.text(paste0(nmsi, "-C", j), x[j], y = cfi[, j], default.units = "native")
      } else {
        grid.text(nmsi, x = x, y = cfi, default.units = "native")
      }
    }
    grid.yaxis(at = c(ceiling(yscale[1] * 100)/100, floor(yscale[2] * 100)/100))
    grid.rect(gp = gpar(fill = "transparent"))

    upViewport(2)
  }
  
  return(panelfun)
}
class(node_profileplot) <- "grapcon_generator"


## region plot visualization function
node_regionplot <- function(mobobj, names = FALSE, abbreviate = TRUE, type = c("mode", "median", "mean"),
  ref = NULL, ylim = NULL, off = 0.1, col_fun = gray.colors,
  uo_show = TRUE, uo_col = "red", uo_lty = 2, uo_lwd = 1.25, ylines = 2)    
{
  ## check input
  stopifnot(!is.null(mobobj))
  stopifnot(off >= 0)
  type <- match.arg(type)

  ## function to extract absolute item threshold parameters from model objects in terminal nodes
  threshparlst <- function (node) threshpar(node, ref = ref, type = type, alias = TRUE, relative = FALSE, cumulative = FALSE, vcov = FALSE)

  ## setup threshold parameters
  node <- nodeids(mobobj, terminal = TRUE)
  delta_lst <- apply_to_models(mobobj, node, FUN = threshparlst, ref = ref, type = type)

  ## setup plotting parameters
  m <- max(sapply(delta_lst, length))
  xi <- 0:m + c(0:(m - 1), m - 1) * off
  xlim <- c(xi[1], xi[m + 1])

  ## setup axis range
  if (is.null(ylim)) ylim <- extendrange(unlist(delta_lst, use.names = FALSE), f = 0.25)

  ## labeling
  if (isTRUE(names)) {
    names <- lapply(delta_lst, names)
  } else if (is.character(names)) {
    names <- split(rep(names, length(node)), f = rep(1:length(node), each = length(names)))
  } else {
    ncf <- lapply(delta_lst, NROW)
    names <- lapply(ncf, function(m) {
      lab <- rep("", m)
      lab[c(1, m)] <- c(1, m)
      pr <- pretty(1:m, n = 4)
      pr <- pr[pr > 1 & pr < m]
      lab[pr] <- pr    
      lab
    })
    abbreviate <- FALSE
  }
  
  ## abbreviation
  if (is.logical(abbreviate)) {
    nlab <- max(unlist(lapply(names, function (j) nchar(j))))
    abbreviate <- if (abbreviate) as.numeric(cut(nlab, c(-Inf, 1.5, 4.5, 7.5, Inf))) else nlab
  }
  names <- lapply(names, function (j) abbreviate(j, abbreviate))

  ## label for extraction
  names(names) <- names(delta_lst) <- node

  ## finally, panel function, just select and paint.
  panelfun <- function(node) {
    ## select node id, coefficients, identified items, x-position vector and label
    id <- as.character(id_node(node))
    delta_unsorted <- delta_lst[[id]]
    lab <- paste("node", id, sep = "")
    namesi <- names[[id]]
    
    ## compute sorted absolute item threshold parameters (code borrowed from psychotools::regionplot())
    delta_sorted <- delta_unsorted
    us <- sapply(delta_unsorted, is.unsorted)
    if (any(us)) {
      usj <- which(us)
      for (j in usj) {
        tpj <- delta_unsorted[[j]]
        nj <- length(tpj)
        
        ## check if there is a point with a biggest parameter, if yes, take mean
        for (i in 1:nj) {
          if (all(tpj[i] > tpj[(i+1):nj])) {
            tpj[i] <- mean(tpj[i:nj])
            tpj <- tpj[-(i+1:nj)]
            break
          }
        }
        
        ## recursive sorting if there is still unorder (e.g. 4, 2, 3, 1)
        while(is.unsorted(tpj)) {
          uo_pos <- which(diff(tpj) < 0)		     # locate unordered parameters, returns position of the first
          tpj[uo_pos] <- (tpj[uo_pos] + tpj[uo_pos + 1]) / 2 # replace first with location of intersection of ccc curves (= (eps1 + eps2)/ 2)
          tpj <- tpj[-(uo_pos + 1)]			     # remove second
        }
        
        delta_sorted[[j]] <- tpj
      }
    }

    ## terminal panel viewport setup
    top.vp <- viewport(layout = grid.layout(nrow = 2, ncol = 1, widths = unit(1, "null"), heights = unit(c(1, 1), c("lines", "null"))),
        		     width = unit(1, "npc"), height = unit(1, "npc") - unit(2, "lines"), name = paste(lab, "_effects", sep = ""))
    pushViewport(top.vp)
    grid.rect(gp = gpar(fill = "white", col = 0), name = paste(lab, "_border", sep = ""))

    ## main title
    pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1, name = paste(lab, "_title_vp", sep = "")))
    grid.text(paste("Node ", id, " (n = ", info_node(node)$nobs, ")", sep = ""), name = paste(lab, "_title", sep = ""))
    upViewport()
    
    ## finally the actual plot
    pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 2, name = lab))
    lab <- paste(lab, "_plot", sep = "")

    ## setup plotting area (3x3)
    wcol <- c(ylines, 1, 1)
    hrow <- c(0.5, 1, 1)

    top.vp <- viewport(layout = grid.layout(nrow = 3, ncol = 3, widths = unit(wcol, c("lines", "null", "lines")), heights = unit(hrow, c("lines", "null", "lines"))), name = paste(lab, "_top_vp", sep = ""))
    bmargin.vp <- viewport(layout.pos.row = 3, layout.pos.col = 2, name = paste(lab, "_bottom-margin_vp", sep = ""))
    lmargin.vp <- viewport(layout.pos.row = 2, layout.pos.col = 1, name = paste(lab, "_left-margin_vp", sep = ""))
    rmargin.vp <- viewport(layout.pos.row = 2, layout.pos.col = 3, name = paste(lab, "_right-margin_vp", sep = ""))
    plot.vp <- viewport(layout.pos.row = 2, layout.pos.col = 2, name = paste(lab, "_vp", sep = ""), xscale = xlim, yscale = ylim)
    pushViewport(top.vp)
    pushViewport(plot.vp)

    ## plot rectangles per item
    for (j in seq_along(delta_sorted)) {
        ncat <- length(delta_sorted[[j]]) + 1
        grid.rect(x = rep.int(xi[j], ncat), y = c(ylim[1], delta_sorted[[j]]), width = rep.int(1, ncat),
        		height = diff.default(c(ylim[1], delta_sorted[[j]], ylim[2])), just = c("left", "bottom"),
        		gp = gpar(fill = col_fun(ncat)), default.units = "native", name = paste(lab, "_item", j, "_rect", sep = ""))
    }

    ## if requested: indicate unordered parameters
    if (uo_show && type == "mode") {
      uo_items <- which(!sapply(mapply(all.equal, delta_sorted, delta_unsorted, check.attributes = FALSE, SIMPLIFY = FALSE, USE.NAMES = FALSE), is.logical))
      for (j in uo_items) {
        uo_pars <- setdiff(delta_unsorted[[j]], delta_sorted[[j]])
        grid.polyline(x = rep(c(xi[j], xi[j] + 1), length(uo_pars)), y = rep(uo_pars, each = 2), default.units = "native", id = rep(1:length(uo_pars), each = 2),
        	      name = paste(lab, "item", j, "_uolines", sep = ""), gp = gpar(col = uo_col, lwd = uo_lwd, lty = uo_lty))
      }
    }

    ## add box and axis
    grid.rect(name = paste(lab, "_plot-box", sep = ""))
    grid.xaxis(at = (xi[-(m+1)] + 0.5), label = namesi, main = TRUE, name = paste(lab, "_xaxis-bottom", sep = ""))
    grid.yaxis(main = TRUE, name = paste(lab, "_yaxis-left", sep = ""))
    upViewport()
    
    ## add descriptions
    pushViewport(lmargin.vp)
    upViewport(2)
    
    ## go back to uper vp
    upViewport(2)
  }

  ## return
  return(panelfun)
}
class(node_regionplot) <- "grapcon_generator"


## bradley-terry plot visualization function
node_btplot <- function(mobobj, id = TRUE, worth = TRUE, names = TRUE,
  abbreviate = TRUE, index = TRUE, ref = TRUE,col = "black", refcol = "lightgray",
  cex = 0.5, pch = 19, xscale = NULL, yscale = NULL, ylines = 1.5)
{
    ## node ids
    node <- nodeids(mobobj, terminal = FALSE)

    ## get all coefficients 
    cf <- apply_to_models(mobobj, node, FUN = function(z)        
      if(worth) worth(z) else coef(z, all = FALSE, ref = TRUE))
    cf <- do.call("rbind", cf)
    rownames(cf) <- node

    ## get one full model
    mod <- apply_to_models(mobobj, node = 1L, FUN = NULL)

    if(!worth) {
      if(is.character(ref) | is.numeric(ref)) {
        reflab <- ref
        ref <- TRUE
      } else {
        reflab <- mod$ref
      }
      if(is.character(reflab)) reflab <- match(reflab, mod$labels)
      cf <- cf - cf[,reflab]
    }

    ## reference
    if(worth) {
      cf_ref <- 1/ncol(cf)
    } else {
      cf_ref <- 0
    }

    ## labeling
    if(is.character(names)) {
      colnames(cf) <- names
      names <- TRUE
    }

    ## abbreviation
    if(is.logical(abbreviate)) {
      nlab <- max(nchar(colnames(cf)))
      abbreviate <- if(abbreviate) as.numeric(cut(nlab, c(-Inf, 1.5, 4.5, 7.5, Inf))) else nlab
    }
    colnames(cf) <- abbreviate(colnames(cf), abbreviate)
    
    if(index) {
      x <- 1:NCOL(cf)
      if(is.null(xscale)) xscale <- range(x) + c(-0.1, 0.1) * diff(range(x))
    } else {
      x <- rep(0, length(cf))
      if(is.null(xscale)) xscale <- c(-1, 1)      
    }
    if(is.null(yscale)) yscale <- range(cf) + c(-0.1, 0.1) * diff(range(cf))
         
    ## panel function for bt plots in nodes
    rval <- function(node) {

      ## node index
      idn <- id_node(node)
    
      ## dependent variable setup
      cfi <- cf[idn, ]

      ## viewport setup
      top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                                   widths = unit(c(ylines, 1, 1), c("lines", "null", "lines")),  
                                   heights = unit(c(1, 1), c("lines", "null"))),
                               width = unit(1, "npc"), 
                               height = unit(1, "npc") - unit(2, "lines"),
                               name = paste("node_btplot", idn, sep = ""))
      pushViewport(top_vp)
      grid.rect(gp = gpar(fill = "white", col = 0))

      ## main title
      top <- viewport(layout.pos.col = 2, layout.pos.row = 1)
      pushViewport(top)
      mainlab <- paste(ifelse(id, paste("Node", idn, "(n = "), ""),
        	       info_node(node)$nobs, ifelse(id, ")", ""), sep = "")
      grid.text(mainlab)
      popViewport()

      ## actual plot  
      plot_vpi <- viewport(layout.pos.col = 2, layout.pos.row = 2,
        xscale = xscale, yscale = yscale,
        name = paste("node_btplot", idn, "plot", sep = ""))
      pushViewport(plot_vpi)

      grid.lines(xscale, c(cf_ref, cf_ref), gp = gpar(col = refcol), default.units = "native")
      if(index) {
        grid.lines(x, cfi, gp = gpar(col = col, lty = 2), default.units = "native")
        grid.points(x, cfi, gp = gpar(col = col, cex = cex), pch = pch, default.units = "native")
        grid.xaxis(at = x, label = if(names) names(cfi) else x)
      } else {  	
        if(names) grid.text(names(cfi), x = x, y = cfi, default.units = "native")
          else grid.points(x, cfi, gp = gpar(col = col, cex = cex), pch = pch, default.units = "native")
      }
      grid.yaxis(at = c(ceiling(yscale[1] * 100)/100, floor(yscale[2] * 100)/100))
      grid.rect(gp = gpar(fill = "transparent"))

      upViewport(2)
    }
	    
    return(rval)
}
class(node_btplot) <- "grapcon_generator"


## MPT plot visualization function
node_mptplot <- function(mobobj, id = TRUE,
  names = TRUE, abbreviate = TRUE, index = TRUE, ref = TRUE,
  col = "black", linecol = "lightgray", cex = 0.5, pch = 19, xscale = NULL,
  yscale = c(0, 1), ylines = 1.5)
{
    ## node ids
    node <- nodeids(mobobj, terminal = FALSE)

    ## get all coefficients
    cf <- apply_to_models(mobobj, node, coef)
    cf <- do.call("rbind", cf)
    rownames(cf) <- node

    ## reference
    if(ref)
      cf_ref <- 1/2

    ## labeling
    if(is.character(names)) {
      colnames(cf) <- names
      names <- TRUE
    }
    if(is.character(index)) {
      cf <- cf[, index]  # reorder labels
      index <- TRUE
    }

    ## abbreviation
    if(is.logical(abbreviate)) {
      nlab <- max(nchar(colnames(cf)))
      abbreviate <- if(abbreviate)
        as.numeric(cut(nlab, c(-Inf, 1.5, 4.5, 7.5, Inf))) else nlab
    }
    colnames(cf) <- abbreviate(colnames(cf), abbreviate)

    if(index) {
      x <- seq_len(NCOL(cf))
      if(is.null(xscale)) {
        if (length(x) > 1) xscale <- range(x) + c(-0.1, 0.1) * diff(range(x))
        else xscale <- c(1 - 0.1, 1 + 0.1)
      }
    } else {
      x <- rep(0, NCOL(cf))
      if(is.null(xscale)) xscale <- c(-1, 1)
    }
    if(is.null(yscale)) {
      if (length(cf) > 1) yscale <- range(cf) + c(-0.1, 0.1) * diff(range(cf))
      else yscale <- c(cf - 0.1, cf + 0.1)
    }

    ## panel function for mpt plots in nodes
    rval <- function(node) {

      ## node index
      idn <- id_node(node)

      ## dependent variable setup
      cfi <- setNames(cf[idn, ], colnames(cf))

      ## viewport setup
      top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                         widths = unit(c(ylines, 1, 1),
                                       c("lines", "null", "lines")),
                         heights = unit(c(1, 1), c("lines", "null"))),
                         width = unit(1, "npc"),
                         height = unit(1, "npc") - unit(2, "lines"),
                         name = paste("node_mptplot", idn, sep = ""))
      pushViewport(top_vp)
      grid.rect(gp = gpar(fill = "white", col = 0))

      ## main title
      top <- viewport(layout.pos.col = 2, layout.pos.row = 1)
      pushViewport(top)
      mainlab <- paste(ifelse(id, paste("Node", idn, "(n = "), ""),
                       info_node(node)$nobs, ifelse(id, ")", ""), sep = "")
      grid.text(mainlab)
      popViewport()

      ## actual plot  
      plot_vpi <- viewport(layout.pos.col = 2, layout.pos.row = 2,
        xscale = xscale, yscale = yscale,
        name = paste("node_mptplot", idn, "plot", sep = ""))
      pushViewport(plot_vpi)

      grid.lines(xscale, c(cf_ref, cf_ref), gp = gpar(col = linecol),
                 default.units = "native")
      if(index) {
        grid.lines(x, cfi, gp = gpar(col = col, lty = 2),
                   default.units = "native")
        grid.points(x, cfi, gp = gpar(col = col, cex = cex), pch = pch,
                    default.units = "native")
        grid.xaxis(at = x, label = if(names) names(cfi) else x)
      } else {
        if(names) grid.text(names(cfi), x = x, y = cfi,
                            default.units = "native")
          else grid.points(x, cfi, gp = gpar(col = col, cex = cex),
                           pch = pch, default.units = "native")
      }
      grid.yaxis(at = c(ceiling(yscale[1] * 100)/100,
                        floor(yscale[2] * 100)/100))
      grid.rect(gp = gpar(fill = "transparent"))

      upViewport(2)
    }

    return(rval)
}
class(node_mptplot) <- "grapcon_generator"

