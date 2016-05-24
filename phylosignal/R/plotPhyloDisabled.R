# This is a slightly modified version of the function plot.phylo
# written by Emmanuel Paradis (Copyright 2002-2014 Emmanuel Paradis)
# and distributed in the ape package (v.3.1-4)
# under the GNU GENERAL PUBLIC LICENSE v.2
.plotPhyloDisabled <- function (x, type = "phylogram", use.edge.length = TRUE, node.pos = NULL, 
                               show.tip.label = TRUE, show.node.label = FALSE, edge.color = "black", 
                               edge.width = 1, edge.lty = 1, font = 3, cex = par("cex"), 
                               adj = NULL, srt = 0, no.margin = FALSE, root.edge = FALSE, 
                               label.offset = 0, underscore = FALSE, x.lim = NULL, y.lim = NULL, 
                               direction = "rightwards", lab4ut = NULL, tip.color = "black", 
                               plot = TRUE, rotate.tree = 0, open.angle = 0, node.depth = 1, 
                               ...) 
{
  Ntip <- length(x$tip.label)
  if (Ntip < 2) {
    warning("found less than 2 tips in the tree")
    return(NULL)
  }
  if (any(tabulate(x$edge[, 1]) == 1)) 
    stop("there are single (non-splitting) nodes in your tree; you may need to use collapse.singles()")
    
  Nedge <- dim(x$edge)[1]
  Nnode <- x$Nnode
  if (any(x$edge < 1) || any(x$edge > Ntip + Nnode)) 
    stop("tree badly conformed; cannot plot. Check the edge matrix.")
  ROOT <- Ntip + 1
  type <- match.arg(type, c("phylogram", "cladogram", "fan", 
                            "unrooted", "radial"))
  direction <- match.arg(direction, c("rightwards", "leftwards", 
                                      "upwards", "downwards"))
  if (is.null(x$edge.length)) 
    use.edge.length <- FALSE
  if (type %in% c("unrooted", "radial") || !use.edge.length || 
        is.null(x$root.edge) || !x$root.edge) 
    root.edge <- FALSE
  if (type == "fan" && root.edge) {
    warning("drawing root edge with type = 'fan' is not yet supported")
    root.edge <- FALSE
  }
  phyloORclado <- type %in% c("phylogram", "cladogram")
  horizontal <- direction %in% c("rightwards", "leftwards")
  xe <- x$edge
  if (phyloORclado) {
    phyOrder <- attr(x, "order")
    if (is.null(phyOrder) || phyOrder != "cladewise") {
      x <- reorder(x)
      if (!identical(x$edge, xe)) {
        ereorder <- match(x$edge[, 2], xe[, 2])
        if (length(edge.color) > 1) {
          edge.color <- rep(edge.color, length.out = Nedge)
          edge.color <- edge.color[ereorder]
        }
        if (length(edge.width) > 1) {
          edge.width <- rep(edge.width, length.out = Nedge)
          edge.width <- edge.width[ereorder]
        }
        if (length(edge.lty) > 1) {
          edge.lty <- rep(edge.lty, length.out = Nedge)
          edge.lty <- edge.lty[ereorder]
        }
      }
    }
    yy <- numeric(Ntip + Nnode)
    TIPS <- x$edge[x$edge[, 2] <= Ntip, 2]
    yy[TIPS] <- 1:Ntip
  }
  z <- reorder(x, order = "postorder")
  if (phyloORclado) {
    if (is.null(node.pos)) 
      node.pos <- if (type == "cladogram" && !use.edge.length) 
        2
    else 1
    if (node.pos == 1) 
      yy <- node.height(x, clado.style = FALSE)
    else {
      xx <- node.depth(x) - 1
      yy <- node.height(x, clado.style = TRUE)
    }
    if (!use.edge.length) {
      if (node.pos != 2) 
        xx <- node.depth(x) - 1
      xx <- max(xx) - xx
    }
    else {
      xx <- node.depth.edgelength(x)
    }
  }
  else {
    twopi <- 2 * pi
    rotate.tree <- twopi * rotate.tree/360
    if (type != "unrooted") {
      TIPS <- x$edge[which(x$edge[, 2] <= Ntip), 2]
      xx <- seq(0, twopi * (1 - 1/Ntip) - twopi * open.angle/360, 
                length.out = Ntip)
      theta <- double(Ntip)
      theta[TIPS] <- xx
      theta <- c(theta, numeric(Nnode))
    }
    switch(type, fan = {
      theta <- node.height(x)
      if (use.edge.length) {
        r <- node.depth.edgelength(x)
      } else {
        r <- node.depth(x)
        r <- 1/r
      }
      theta <- theta + rotate.tree
      xx <- r * cos(theta)
      yy <- r * sin(theta)
    }, unrooted = {
      nb.sp <- node.depth(x)
      XY <- if (use.edge.length) unrooted.xy(Ntip, Nnode, 
                                             z$edge, z$edge.length, nb.sp, rotate.tree) else unrooted.xy(Ntip, 
                                                                                                         Nnode, z$edge, rep(1, Nedge), nb.sp, rotate.tree)
      xx <- XY$M[, 1] - min(XY$M[, 1])
      yy <- XY$M[, 2] - min(XY$M[, 2])
    }, radial = {
      r <- node.depth(x)
      r[r == 1] <- 0
      r <- 1 - r/Ntip
      theta <- node.height(x) + rotate.tree
      xx <- r * cos(theta)
      yy <- r * sin(theta)
    })
  }
  if (phyloORclado) {
    if (!horizontal) {
      tmp <- yy
      yy <- xx
      xx <- tmp - min(tmp) + 1
    }
    if (root.edge) {
      if (direction == "rightwards") 
        xx <- xx + x$root.edge
      if (direction == "upwards") 
        yy <- yy + x$root.edge
    }
  }
  if (no.margin) 
    par(mai = rep(0, 4))
  if (show.tip.label) 
    nchar.tip.label <- nchar(x$tip.label)
  max.yy <- max(yy)
  if (is.null(x.lim)) {
    if (phyloORclado) {
      if (horizontal) {
        x.lim <- c(0, NA)
        pin1 <- par("pin")[1]
        strWi <- strwidth(x$tip.label, "inches", cex = cex)
        xx.tips <- xx[1:Ntip] * 1.04
        alp <- try(uniroot(function(a) max(a * xx.tips + 
                                             strWi) - pin1, c(0, 1e+06))$root, silent = TRUE)
        if (is.character(alp)) 
          tmp <- max(xx.tips) * 1.5
        else {
          tmp <- if (show.tip.label) 
            max(xx.tips + strWi/alp)
          else max(xx.tips)
        }
        if (show.tip.label) 
          tmp <- tmp + label.offset
        x.lim[2] <- tmp
      }
      else x.lim <- c(1, Ntip)
    }
    else switch(type, fan = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy * 
                        cex)
        x.lim <- range(xx) + c(-offset, offset)
      } else x.lim <- range(xx)
    }, unrooted = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy * 
                        cex)
        x.lim <- c(0 - offset, max(xx) + offset)
      } else x.lim <- c(0, max(xx))
    }, radial = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.03 * cex)
        x.lim <- c(-1 - offset, 1 + offset)
      } else x.lim <- c(-1, 1)
    })
  }
  else if (length(x.lim) == 1) {
    x.lim <- c(0, x.lim)
    if (phyloORclado && !horizontal) 
      x.lim[1] <- 1
    if (type %in% c("fan", "unrooted") && show.tip.label) 
      x.lim[1] <- -max(nchar.tip.label * 0.018 * max.yy * 
                         cex)
    if (type == "radial") 
      x.lim[1] <- if (show.tip.label) 
        -1 - max(nchar.tip.label * 0.03 * cex)
    else -1
  }
  if (phyloORclado && direction == "leftwards") 
    xx <- x.lim[2] - xx
  if (is.null(y.lim)) {
    if (phyloORclado) {
      if (horizontal) 
        y.lim <- c(1, Ntip)
      else {
        y.lim <- c(0, NA)
        pin2 <- par("pin")[2]
        strWi <- strwidth(x$tip.label, "inches", cex = cex)
        yy.tips <- yy[1:Ntip] * 1.04
        alp <- try(uniroot(function(a) max(a * yy.tips + 
                                             strWi) - pin2, c(0, 1e+06))$root, silent = TRUE)
        if (is.character(alp)) 
          tmp <- max(yy.tips) * 1.5
        else {
          tmp <- if (show.tip.label) 
            max(yy.tips + strWi/alp)
          else max(yy.tips)
        }
        if (show.tip.label) 
          tmp <- tmp + label.offset
        y.lim[2] <- tmp
      }
    }
    else switch(type, fan = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy * 
                        cex)
        y.lim <- c(min(yy) - offset, max.yy + offset)
      } else y.lim <- c(min(yy), max.yy)
    }, unrooted = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy * 
                        cex)
        y.lim <- c(0 - offset, max.yy + offset)
      } else y.lim <- c(0, max.yy)
    }, radial = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.03 * cex)
        y.lim <- c(-1 - offset, 1 + offset)
      } else y.lim <- c(-1, 1)
    })
  }
  else if (length(y.lim) == 1) {
    y.lim <- c(0, y.lim)
    if (phyloORclado && horizontal) 
      y.lim[1] <- 1
    if (type %in% c("fan", "unrooted") && show.tip.label) 
      y.lim[1] <- -max(nchar.tip.label * 0.018 * max.yy * 
                         cex)
    if (type == "radial") 
      y.lim[1] <- if (show.tip.label) 
        -1 - max(nchar.tip.label * 0.018 * max.yy * cex)
    else -1
  }
  if (phyloORclado && direction == "downwards") 
    yy <- y.lim[2] - yy
  if (phyloORclado && root.edge) {
    if (direction == "leftwards") 
      x.lim[2] <- x.lim[2] + x$root.edge
    if (direction == "downwards") 
      y.lim[2] <- y.lim[2] + x$root.edge
  }
  asp <- if (type %in% c("fan", "radial", "unrooted")) 
    1
  else NA
  
  L <- list(type = type, use.edge.length = use.edge.length, 
            node.pos = node.pos, node.depth = node.depth, show.tip.label = show.tip.label, 
            show.node.label = show.node.label, font = font, cex = cex, 
            adj = adj, srt = srt, no.margin = no.margin, label.offset = label.offset, 
            x.lim = x.lim, y.lim = y.lim, direction = direction, 
            tip.color = tip.color, Ntip = Ntip, Nnode = Nnode)
  assign("last_plot.phylo", c(L, list(edge = xe, xx = xx, yy = yy)), 
         envir = .PlotPhyloEnv)
  return(L)
}