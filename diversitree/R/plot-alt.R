## New "alternative" plotting interface for diversitree.
## TODO: node labels not done.

## The argument list here is exactly the same as that for
## ape:::plot.phylo() - additional arguments will go in after ...?
plot2.phylo <- function(x, type="phylogram", use.edge.length=TRUE,
                        node.pos=NULL, show.tip.label=TRUE,
                        show.node.label=FALSE, edge.color="black",
                        edge.width=1, edge.lty=1, font=3,
                        cex=par("cex"), adj=NULL, srt=0,
                        no.margin=FALSE, root.edge=FALSE,
                        label.offset=0, underscore=FALSE,
                        x.lim=NULL, y.lim=NULL,
                        direction="rightwards", lab4ut="horizontal",
                        tip.color="black", ...,
                        n.taxa=NULL,
                        clade.lwd=1,
                        clade.color=NULL, clade.fill=NULL, pad=0) {
  n.tip <- length(x$tip.label)
  ## Difference: ape has warning and returns NULL - I have error.
  if ( n.tip < 2 )
    stop("Cannot plot tree with < 2 tips")
  if ( any(tabulate(x$edge[, 1]) == 1) )
    stop("there are single (non-splitting) nodes in your tree;\n",
         "\tYou may need to use collapse.singles()")

  ## TODO: cladogram, fan, unrooted are not implemented.
  type <- match.arg(type, c("phylogram", "fan"))
  ## TODO: other directions are not implemented.
  direction <- match.arg(direction, "rightwards")

  if ( is.null(x$edge.length) )
    use.edge.length <- FALSE

  if ( no.margin )
    par(mar=rep(0, 4))

  ## Difference: I think this is handled differently by ape
  ## Repeat arguments appropriately - all called functions may rely on
  ## these being appropriately vectorised.
  tip.color <- rep(tip.color, length.out=n.tip)
  edge.color <- rep(edge.color, length.out=nrow(x$edge))
  edge.width <- rep(edge.width, length.out=nrow(x$edge))
  edge.lty <- rep(edge.lty, length.out=nrow(x$edge))

  ## TODO: the clade borders need repeating here?

  ## Difference: No ape analogue here
  if ( !is.null(n.taxa) ) {
    if ( !is.null(names(n.taxa)) ) {
      n.taxa <- n.taxa[x$tip.label]
      n.taxa[is.na(n.taxa)] <- 1
    } else if ( length(n.taxa) != n.tip ) {
      stop("n.taxa must be same length as x$tip.label")
    }
  }
  x$n.taxa <- n.taxa
  ## This is used to compute the y limits, and turns up over and
  ## over.
  x$n.spp <- length(x$tip.label) + sum(n.taxa - 1)

  ## This is not perfect, but should cover most cases.
  if ( !is.null(x$n.taxa) ) {
    i <- match(seq_len(n.tip), x$edge[,2])
    rep.clade <- function(x, edge) {
      if (is.null(x))
        x <- edge
      else if (length(x) == 1)
        x <- rep(x, n.tip)
      else if (length(x) != n.tip)
        stop("Not yet handled")
      x
    }

    clade.color <- rep.clade(clade.color, edge.color[i])
    clade.fill  <- rep.clade(clade.fill,  edge.color[i])
    clade.lwd   <- rep.clade(clade.lwd,   edge.width[i])

    ## Supress drawing of the terminals that have clades.
    edge.color[match(which(x$n.taxa > 1), x$edge[,2])] <- NA
  }  

  ## All the plotting functions are fairly similar, differing mostly
  ## in the exact name of the functions used, but in about the same
  ## order.
  pp.coords <- get(paste("pp.coords", type, sep="."))
  pp.lim <- get(paste("pp.lim", type, sep="."))
  pp.segments <- get(paste("pp.segments", type, sep="."))
  pp.tiplabel <- get(paste("pp.tiplabel", type, sep="."))
  pp.clades <- get(paste("pp.clades", type, sep="."))

  asp <- if ( type == "fan" ) 1 else NA

  if ( is.null(node.pos) )
    node.pos <- 1
  else if ( node.pos != 1 )
    stop("node.pos != 1 not yet implemented")

  if ( is.null(adj) )
    adj <- 0
  
  xy <- pp.node.coords(x)
  xy.seg <- pp.coords(x, xy)
 
  lims <- pp.lim(x, xy, x.lim, y.lim, cex,
                 show.tip.label, label.offset + pad)
  plot(NA, type="n", xlim=lims$xlim, ylim=lims$ylim, xlab="",
       ylab="", xaxt="n", yaxt="n", bty="n", asp=asp, ...)
  pp.segments(x, xy.seg, edge.color, edge.width, edge.lty)

  if ( !is.null(x$n.taxa) && any(x$n.taxa > 1) )
    pp.clades(x, xy, xy.seg, clade.color, clade.fill, clade.lwd)
  if ( show.tip.label ) {
    if ( !underscore )
      x$tip.label <- gsub("_", " ", x$tip.label)
    pp.tiplabel(x, xy, label.offset, adj, cex, tip.color, font)
  }

  if ( type == "fan" )
    xy <- pp.coords.fix.xy(x, xy)

  ## This has all the arguments below for compatibility with ape.
  ret <- list(type=type, use.edge.length=use.edge.length,
              node.pos=node.pos, show.tip.label=show.tip.label,
              show.node.label=show.node.label, font=font,
              cex=cex, adj=adj, srt=srt, no.margin=no.margin,
              label.offset=label.offset,
              x.lim=lims$xlim, y.lim=lims$ylim,
              direction=direction, tip.color=tip.color,
              Ntip=n.tip, Nnode=x$Nnode,
              edge=x$edge,
              xx=xy$xx,
              yy=xy$yy,
              ## Extra below here
              n.taxa=x$n.taxa, n.spp=x$n.spp, xy=xy, xy.seg=xy.seg)
  assign("last_plot.phylo", ret, envir=.PlotPhyloEnv)
  invisible(ret)
}

## Returns x/y coordinates of the nodes in a phylogeny.
## This is derived from ape's nodeDepth
pp.node.coords <- function(phy) {
  phy.p <- reorder(phy, "pruningwise")
  if ( is.null(phy.p$edge.length) )
    xx <- node.depth(phy.p)
  else
    xx <- node.depth.edgelength(phy.p)
  yy <- node.height(phy.p, phy$n.taxa)
  data.frame(xx=xx, yy=yy)
}

## This returns a matrix with x/y coordinates for the different
## segments of a tree.
pp.coords.phylogram <- function(phy, xy) {
  phy.p <- reorder(phy, "pruningwise")
  phy.c <- reorder(phy, "cladewise")
  
  edge <- phy.c$edge
  n.node <- phy$Nnode
  n.tip <- length(phy$tip.label)
  nodes <- (n.tip + 1):(n.tip + n.node)

  ## First, grab the easy ones:
  x0v <- xy$xx[nodes]
  y0v <- y1v <- numeric(n.node)
  x0h <- xy$xx[edge[,1]]
  x1h <- xy$xx[edge[,2]]
  y0h <- xy$yy[edge[,2]]

  for ( i in nodes ) {
    tmp <- range(xy$yy[edge[edge[,1] == i, 2]])
    y0v[i - n.tip] <- tmp[1]
    y1v[i - n.tip] <- tmp[2]
  }

  h <- data.frame(x0=x0h, y0=y0h, x1=x1h, y1=y0h, horiz=TRUE,
                  idx=seq_len(nrow(edge)) + n.tip + n.node,
                  parent=edge[,1])
  v <- data.frame(x0=x0v, y0=y0v, x1=x0v, y1=y1v, horiz=FALSE,
                  idx=seq_len(n.node) + n.tip,
                  parent=h$idx[match(nodes, edge[,2])])
  rbind(v, h)
}

pp.coords.fan <- function(phy, xy) {
  xy.seg <- pp.coords.phylogram(phy, xy)
  
  n <- length(phy$tip.label) + sum(phy$n.taxa - 1)
  xy.seg$theta0 <- xy.seg$y0 / (n + 1) * 2 * pi
  xy.seg$theta1 <- xy.seg$y1 / (n + 1) * 2 * pi
  xy.seg$r0 <- xy.seg$x0
  xy.seg$r1 <- xy.seg$x1
  xy.seg
}

pp.coords.fix.xy <- function(phy, xy) {
  n <- length(phy$tip.label) + sum(phy$n.taxa - 1)  
  xy$theta <- xy$yy / (n + 1) * 2 * pi
  xy$r <- xy$xx

  xy$xx <- with(xy, r * cos(theta))
  xy$yy <- with(xy, r * sin(theta))

  xy
}

## Compute limits for the plot.
pp.lim.phylogram <- function(x, xy, xlim, ylim, cex,
                             show.tip.label, label.offset) {
  n.tip <- length(x$tip.label)
  if ( is.null(xlim) ) {
    xlim <- c(0, NA)
    xlab <- if ( show.tip.label ) 
      nchar(x$tip.label) * 0.018 * max(xy$xx) * cex else 0
    xlim[2] <- max(xy$xx[1:n.tip] + xlab)
    xlim[2] <- xlim[2] + label.offset
  }
  if ( is.null(ylim) ) {
    if ( is.null(x$n.taxa) )
      ylim <- c(1, n.tip)
    else
      ylim <- c(1, x$n.spp)
  }

  list(xlim=xlim, ylim=ylim)  
}

pp.lim.fan <- function(x, xy, xlim, ylim, cex,
                          show.tip.label, label.offset) {
  xmax <- max(xy$xx) + label.offset

  ## TODO: Fix this by computing what r + strwidth(.) + label.offset
  ## is, then transforming into cartesian coordinates?  But then have
  ## to backtrack a bunch as usr coordinates not set up yet.
  if ( show.tip.label ) {
    pin <- min(par("pin"))
    lab <- max(strwidth(x$tip.label, "inches", cex)) * 2
    if ( lab > pin )
      stop("Margins too large - reduce cex and/or margins")
    xmax <- xmax / (1 - lab / pin)
  }
  
  if ( is.null(xlim) )
    xlim <- c(-xmax, xmax)
  if ( is.null(ylim) )
    ylim <- c(-xmax, xmax)
  list(xlim=xlim, ylim=ylim)
}

## Simplify adding the line segments to the plot.
pp.segments.phylogram <- function(phy, xy.seg, col, lwd, lty) {
  i <- match(seq_len(phy$Nnode) + length(phy$tip.label), phy$edge[,2])
  j <- c(i, seq_along(col))
  j[is.na(j)] <- 1 # TODO: This should be done better.
  segments(xy.seg$x0, xy.seg$y0, xy.seg$x1, xy.seg$y1,
           col=col[j], lwd=lwd[j], lty=lty[j])
}

## TODO: There is no way of modifying the np here.
pp.segments.fan <- function(phy, xy.seg, col, lwd, lty, np=1000) {
  if ( length(col) == 1 ) col <- rep(col, nrow(phy$edge))
  if ( length(lwd) == 1 ) lwd <- rep(lwd, nrow(phy$edge))
  if ( length(lty) == 1 ) lty <- rep(lty, nrow(phy$edge))
  xy.seg2 <- split(xy.seg, xy.seg$horiz)
  with(xy.seg2[[2]],
       segments(r0*cos(theta0), r0*sin(theta0),
                r1*cos(theta1), r1*sin(theta1),
                lwd=lwd, col=col))

  i <- match(seq_len(phy$Nnode) + length(phy$tip.label), phy$edge[,2])
  if ( any(!is.na(i)) ) # For a two branch tree.
    with(xy.seg2[[1]],
         arcs(theta0, theta1, r1, col[i], lty[i], lwd[i], np))
}  

pp.clades.phylogram <- function(phy, xy, xy.seg, border, fill, lwd) {
  n <- phy$n.spp
  i <- phy$n.taxa > 1
  n.taxa <- phy$n.taxa[i]

  dy <- (n.taxa/2 - .5)
  
  tmp <- xy.seg[xy.seg$horiz,][match(which(i), phy$edge[,2]),] 

  x0 <- tmp$x0
  x1 <- tmp$x1
  ym <- tmp$y0
  y0 <- ym - dy
  y1 <- ym + dy

  graphics::polygon(rbind(x0, x1, x1, x0, NA), 
          rbind(ym, y0, y1, ym, NA),
          border=border[i], col=fill[i], lwd=lwd[i])
}

pp.clades.fan <- function(phy, xy, xy.seg, border, fill, lwd) {
  n <- phy$n.spp
  i <- phy$n.taxa > 1
  n.taxa <- phy$n.taxa[i]
  
  dt <- (n.taxa/2 - .5) / n * 2 * pi

  tmp <- xy.seg[xy.seg$horiz,][match(which(i), phy$edge[,2]),]
  r0 <- tmp$r0
  r1 <- tmp$r1
  t0 <- tmp$theta0 - dt
  t1 <- tmp$theta0 + dt

  sectors(t0, t1, r0, r1, border=border[i], col=fill[i], lwd=lwd)
}

pp.tiplabel.phylogram <- function(phy, xy, label.offset, adj, cex,
                                  col, font) {
  if ( is.null(adj) ) adj <- 0
  n.tip <- length(phy$tip.label)
  wmax <- max(strwidth(phy$tip.label, cex=cex))
  text(xy$xx[1:n.tip] + label.offset + wmax * 1.05 * adj,
       xy$yy[1:n.tip], phy$tip.label, adj=adj, cex=cex, col=col,
       font=font)
}

pp.tiplabel.fan <- function(phy, xy, label.offset, adj, cex,
                            col, font) {
  ## TODO: adj is ignored - should I warn?
  n.tip <- length(phy$tip.label)
  n <- n.tip + sum(phy$n.taxa - 1)
  r     <- xy$xx[1:n.tip] + label.offset
  theta <- xy$yy[1:n.tip] / (n + 1) * 2 * pi
  radial.text(r, theta, phy$tip.label, cex=cex, col=col, font=font)
}


###########################################################################
## Utility functions.  These replace some C functions within ape.
node.depth <- function(x) {
  x <- reorder(x, "pruningwise")  
  n.tip <- length(x$tip.label)
  n.edge <- nrow(x$edge)
  n.node <- x$Nnode
  edge1 <- x$edge[,1]
  edge2 <- x$edge[,2]

  xx <- numeric(n.tip + n.node)
  xx[1:n.tip] <- 1
  for ( i in 1:n.edge)
    xx[edge1[i]] <- xx[edge1[i]] + xx[edge2[i]]
  xx <- xx - 1
  max(xx) - xx
}

node.depth.edgelength <- function(x) {
  x <- reorder(x, "pruningwise")
  n.tip <- length(x$tip.label)
  n.edge <- nrow(x$edge)
  n.node <- x$Nnode
  edge1 <- x$edge[,1]
  edge2 <- x$edge[,2]
  edge.length <- x$edge.length

  xx <- numeric(n.tip + n.node)
  for ( i in n.edge:1 )
    xx[edge2[i]] <- xx[edge1[i]] + edge.length[i]
  xx
}

node.height <- function(x, n.taxa) {
  n.tip <- length(x$tip.label)
  n.edge <- nrow(x$edge)
  n.node <- x$Nnode

  if (!is.null(attr(x, "order")) && attr(x, "order") == "pruningwise")
    x <- reorder(x)

  xx <- reorder(x, order="pruningwise")
  edge1 <- xx$edge[,1]
  edge2 <- xx$edge[,2]

  yy <- numeric(n.tip + n.node)
  tips <- x$edge[x$edge[, 2] <= n.tip, 2]

  if ( !is.null(n.taxa) ) {
    ## First, ensure that all tips (even those that represent singles)
    ## have count information:
    i <- match(names(n.taxa), x$tip.label)
    n.taxa2 <- rep(1, n.tip)
    n.taxa2[i] <- n.taxa

    ## Next, organise things into the plotting order
    n.taxa.plt <- n.taxa2[tips]
    yy[seq_len(n.tip)] <-
      (cumsum(n.taxa.plt) - (n.taxa.plt-1)/2)[order(tips)]
  } else {
    yy[tips] <- 1:n.tip
  }

  S <- n <- 0
  for ( i in 1:n.edge ) {
    S <- S + yy[edge2[i]]
    n <- n + 1
    if ( i == n.edge || edge1[i+1] != edge1[i] ) {
      yy[edge1[i]] <- S/n
      S <- n <- 0
    }
  }
  yy
}


## These things are not properly checked:
pp.nodepoints.phylogram <- function (x, xy, pch=19, ...) {
  points(as.data.frame(xy)[-(1:length(x$tip.label)), ], pch=pch, 
         ...)
}
pp.nodelabel.phylogram <- function(x, xy, label.offset, ...) {
  root <- length(x$tip.label) + 1
  text(xy$xx[root:length(xy$xx)] + label.offset, xy$yy[root:length(xy$yy)], 
       x$node.label, ...)
}
pp.tippoints.phylogram <- function (x, xy, pch=19, ...) {
  with(as.data.frame(xy)[1:length(x$tip.label), ],
       points(xx + 0.5, yy, pch=pch, ...))
}
