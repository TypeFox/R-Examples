## A history object is a list with elements:
##   tip.state
##   node.state
##   discrete - is this a discrete or continuous valued history
##   states - possible values (d) or range (c) of history values
## and possibly
##   history
##   phy
## Other elements could be added at will.
## The format of the 'history' element is a list with one element per
## edge.  Each element contains a matrix where the first column is
## time, measured from the start of the branch (t=0) and the second
## column is the state at that time.
make.history <- function(phy, tip.state, node.state, history,
                         discrete=TRUE, states=NULL, check=TRUE) {
  if ( check ) {
    if ( length(tip.state) != length(phy$tip.label) )
      stop("tip.state of wrong length")
    if ( length(node.state) != phy$Nnode )
      stop("node.state of wrong length")
    if ( length(history) != nrow(phy$edge) )
      stop("history of wrong length")
  }
  structure(list(tip.state=tip.state,
                 node.state=node.state,
                 history=history,
                 discrete=discrete,
                 states=states),
            class="history")
}

history.from.sim.discrete <- function(phy, states) {
  hist <- phy$hist
  if ( is.null(hist) )
    stop("Cannot make a history object")
  ## First, produce a vector of tip and node states in the same order
  ## as the edge matrix.  If no change occured along an edge, then the
  ## base state is the same as the tip state, so assume that this is
  ## correct along a branch
  all.names <- c(phy$tip.label, phy$node.label)[phy$edge[,2]]
  all.states <- c(phy$tip.state, phy$node.state)
  h <- lapply(all.states[all.names],
              function(x) matrix(c(0, x), 1, 2))

  ## A few edges have character changes along them; these are the
  ## elements in 'hist' (h.c: history-changing).  From the information
  ## contained in the phylogeny's history store, we know what time the
  ## event occured at (t), what time the branch started at (x0) and
  ## what time along the branch the change occured at (tc=t-x0).
  h.c <- lapply(split(hist, hist$name2), function(x)
                cbind(c(0, x$tc), c(x$from[1], x$to)))
  h[names(h.c)] <- h.c

  make.history(phy, phy$tip.state, phy$node.state, h, TRUE, states)
}

## Functions involved with plotting histories.
## The only thing 'discrete' specific about this is that the state ->
## colour mapping is assumed to be 1:1.  If the state information is
## easy to map to a colour, this will require only trivial changes to
## make into a full continuous trait mapping.  However, I need this
## for a mapping for (say) summaries over stochastic realisations, so
## we shall see.
plot.history <- function(x, phy, cols=seq_along(states),
                         states=x$states,
                         xlim=NULL, ylim=NULL,
                         show.tip.label=TRUE,
                         show.node.label=FALSE,
                         show.tip.state=TRUE,
                         show.node.state=TRUE,
                         no.margin=FALSE,
                         cex=1, font=3,
                         srt=0, adj=0,
                         label.offset=NA,
                         lwd=1,
                         ...) {
  if ( is.na(label.offset) )
    label.offset <- if ( show.tip.state ) 1 else 0
  discrete <- x$discrete

  xy <- pp.node.coords(phy)
  obj <- plot.history.coords(phy, xy, x)

  if ( no.margin )
    par(mar=rep(0, 4))

  lims <- pp.lim.phylogram(phy, xy, xlim, ylim, cex,
                           show.tip.label, label.offset)
  plot(NA, type="n", xlim=lims$xlim, ylim=lims$ylim, xlab="",
       ylab="", xaxt="n", yaxt="n", bty="n", asp=NA, ...)

  
  ## Improve the colour->tip mapping.  There are a couple of options -
  ## we could have data in {0,1} or {1,2,...,k}
  ## Both may be possibly missing states.
  if ( !discrete && missing(cols) )
    cols <- heat.colors(12)
  if ( (discrete && length(cols) < length(states)) ||
       (!discrete && length(cols) < 2) )
    stop("'cols' too short")

  if ( discrete ) {
    cols.seg <- cols[match(obj$state, states)]
    cols.node <- cols[match(x$node.state, states)]
    cols.tip <- cols[match(x$tip.state, states)]
  } else {
    obj$state.bin <- discretize(obj$state, length(cols), states)
    obj <- history.bin(obj)
    cols.seg <- cols[discretize(obj$state, length(cols), states)]
    cols.node <- cols[discretize(x$node.state, length(cols), states)]
    cols.tip <- cols[discretize(x$tip.state, length(cols), states)]
  }

  with(obj, segments(x0, y0, x1, y1, col=cols.seg, lwd=lwd))

  if ( show.tip.label )
    pp.tiplabel.phylogram(phy, xy, label.offset, col="black",
                          cex=cex, font=font, adj=adj)
  if ( show.node.label )
    pp.nodelabel.phylogram(phy, xy, label.offset,
                           cex=cex, font=font, srt=srt, adj=adj)
  if ( show.node.state )
    pp.nodepoints.phylogram(phy, xy, cex=cex, col=cols.node)
  if ( show.tip.state )
    pp.tippoints.phylogram(phy, xy, cex=cex, col=cols.tip)


  ## For ape compatibility.
  ret <- list(type="cladogram", use.edge.length=TRUE,
              node.pos=NULL, show.tip.label=show.tip.label,
              show.node.label=show.node.label, font=font,
              cex=cex, adj=adj, srt=srt, no.margin=no.margin,
              label.offset=label.offset,
              x.lim=lims$xlim, y.lim=lims$ylim,
              direction="rightwards", tip.color="black",
              Ntip=length(phy$tip.label), Nnode=phy$Nnode,
              edge=phy$edge, xx=xy$xx, yy=xy$yy,
              ## Extra below here
              n.taxa=phy$n.taxa, n.spp=phy$n.spp, xy=xy, xy.seg=obj)
  assign("last_plot.phylo", ret, envir=.PlotPhyloEnv)
  invisible(ret)
}

## This generates the segment matrix, as in plot.phylo.coords, but
## this may include multiple segments per horizontal branch.  State
## information along each branch is also included in the 'state'
## column.
plot.history.coords <- function(phy, xy, hist) {
  s <- hist$node.state
  h <- hist$history

  obj <- pp.coords.phylogram(phy, xy)

  ## First, get the state at the *base* of each branch
  obj$state <- c(s, sapply(h, function(x) x[1,2]))

  ## Check that this is correct - good.
  ##   all.equal(plot.phylo.node.state(reorder(phy, "cladewise"),
  ##                                   subset(obj, horiz)$state),
  ##             as.numeric(s))

  ## Find the branches with more than one state; 'i' is their index
  ## and 'n' is the number of states.
  n <- sapply(h, nrow)
  i <- which(n > 1)
  if ( length(i) > 0 ) {
    n <- n[i]

    ## Construct a new bit of plotting for this branch:
    tmp <- obj[rep(i + phy$Nnode, n),]
    rownames(tmp) <- NULL
    new <- do.call(rbind, h[i])
    tmp$x0 <- tmp$x0 + new[,1]
    j <- unlist(lapply(n, function(i) rep(c(FALSE, TRUE), c(1, i-1))))
    tmp$x1[which(j)-1] <- tmp$x0[j]
    tmp$state <- new[,2]

    ## Add this into the plotting object
    obj[c(i + phy$Nnode, seq_len(sum(n-1)) + nrow(obj)),] <- tmp
    rownames(obj) <- NULL
  }

  obj
}

## It would be nice to offer time-based binning too...
history.bin <- function(obj) {
  f <- function(x) {
    if ( nrow(x) > 1 ) {
      i <- cumsum(c(0, rle(x$state.bin)$lengths))
      z <- x[i[-length(i)] + 1,]
      z$x1 <- x$x1[i[-1]]
      z
    } else {
      x
    }
  }

  horiz <- obj$horiz # R CMD check NOTE avoidance
  tmp <- subset(obj, horiz)
  ret <- rbind(subset(obj, !horiz),
                do.call(rbind, lapply(split(tmp, tmp$idx), f)))
  rownames(ret) <- NULL
  ret
}
