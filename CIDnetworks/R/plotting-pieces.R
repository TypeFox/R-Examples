
# Graphical components common to lots of different subclasses.

###############################################################
# Whole network plots.


# image plot for a network.
image.netplot <- function (edge.list, outcome,
                     extremes=range(outcome,na.rm=TRUE),
                     colvalues=10,
                     col1="white", col2="black",
                     null.color="#DDDDFF",
                     colpalette=colorRampPalette(c(col1,col2)),
                     n.nodes=max(c(edge.list)),
                     lwd=0.1,

                     node.labels=1:n.nodes,
                     label.cex=10/n.nodes,

                     node.order=1:n.nodes,
                     ...) {
  #edge.list=make.edge.list(10); outcome=seq(0, 1, length=nrow(edge.list)); extremes=range(outcome); colvalues=2; col1="white"; col2="black"; colpalette=colorRampPalette(c(col1,col2)); n.nodes=max(c(edge.list)); symmetric=TRUE; null.color="#004400"

  ##  symmetric=TRUE,
  val1 <- unique(edge.list)
  val2 <- unique(rbind(val1, val1[,2:1]))

  ## is there any overlap? If not, symmetric plot.
  symmetric <- (nrow(val2) == 2*nrow(val1))

  plot(c(-1, n.nodes), c(-1, n.nodes), ty="n", axes=FALSE, xlab="Receivers", ylab="Senders", ...)
  rect(0, 0, n.nodes, n.nodes, col=null.color)
  colseq <- colpalette(colvalues)[as.numeric(cut(outcome, breaks=seq(min(extremes), max(extremes), length=colvalues+1), include.lowest=TRUE))]
  colseq[is.na(colseq)] <- "pink"
  f.order <- rep(NA,length(node.order))
  for(ii in 1:length(node.order)){
      f.order[node.order[ii]] <- ii
  }
  node.order <- f.order
  reordered.edge.list <- array(node.order[edge.list], dim(edge.list))

  #rect(reordered.edge.list[,1]-1, reordered.edge.list[,2]-1,
  #     reordered.edge.list[,1], reordered.edge.list[,2], col=colseq, lwd=lwd)
  #if (symmetric) rect(reordered.edge.list[,2]-1, reordered.edge.list[,1]-1,
  #                    reordered.edge.list[,2], reordered.edge.list[,1], col=colseq, lwd=lwd)

  rect(reordered.edge.list[,2]-1, reordered.edge.list[,1]-1,
       reordered.edge.list[,2], reordered.edge.list[,1], col=colseq, lwd=lwd)
  if (symmetric) rect(reordered.edge.list[,1]-1, reordered.edge.list[,2]-1,
                      reordered.edge.list[,1], reordered.edge.list[,2], col=colseq, lwd=lwd)

  text(1:n.nodes - 0.5, rep(-0.5, length(n.nodes)), node.labels[node.order], cex=label.cex)
  text(rep(0, length(n.nodes)), 1:n.nodes - 0.5, node.labels[node.order], cex=label.cex, pos=2)

}


###############################################################
# Latent space plots.

# Projection of points into 2 dimensions for latent space/vector.

project.into.2d <- function(point.matrix) {  #which is n-by-k.
  #point.matrix=matrix(rnorm(4*900), ncol=4)

  if (ncol(point.matrix) > 2) {
    new.2d <- princomp(point.matrix)$scores[,1:2]
  } else if (ncol(point.matrix) == 2) {new.2d <- point.matrix} else
    new.2d <- cbind(point.matrix, point.matrix)/sqrt(2)

  return(new.2d)

}

latent.space.plot <- function (point.matrix, arrowlines=FALSE, labels=1:nrow(point.matrix), ...) {
  #point.matrix=matrix(rnorm(4*900), ncol=4)

  n.nodes <- nrow(point.matrix)
  twod.points <- project.into.2d (point.matrix)
  plot(twod.points, ty="n", ...)
  text(twod.points[,1], twod.points[,2], labels)

  if (arrowlines) segments (rep(0, n.nodes),
                            rep(0, n.nodes),
                            twod.points[,1],
                            twod.points[,2])

}


###############################################################
# Block membership plots.

number.to.vector <- function(number.membership.vector, n.groups=max(number.membership.vector))
  sapply(number.membership.vector, function(kk) {out <- rep(0, n.groups); out[kk] <- 1; out})


color.block.one <- function(color.matrix) {
  #color.matrix = -4:4
  #print(color.matrix)

  red <- 1-abs(color.matrix)/max(6, abs(color.matrix))*(color.matrix>0)
  blue <- 1-abs(color.matrix)/max(6, abs(color.matrix))*(color.matrix<0)
  green <- red*(color.matrix>0) + blue*(color.matrix<0) + 1*(color.matrix==0)

  rgb(red, green, blue)
}

block.membership.plot <- function (membership.share,
                                   block.matrix,
                                   main="Stochastic Block Model Outcome",
                                   node.labels=1:dim(membership.share)[2], ...) {

  #membership.share=matrix(runif(20*2), nrow=2); block.matrix=matrix(2*runif(2*2)-1, nrow=2); main="Stochastic Block Model Outcome"
  #print(block.matrix)

  #what's the minimum number
  n.nodes <- ncol(membership.share)
  n.groups <- nrow(block.matrix)  #
  n.groups.known <- nrow(membership.share)
  min.x <- n.groups + n.nodes + 2
  min.y <- n.groups + 2
  width <- ceiling(sqrt(min.x*min.y))
  #what's the total number of x units rounded? Extend x by 1 until the reserve is enough for the block.
  unit.cap <- (width - (n.nodes %% width)) %% width
  while (unit.cap < n.groups + 2 & width < min.x) {
    width <- width+1
    unit.cap <- (width - (n.nodes %% width) %% width)
  }
  #plot the membership share matrix.

  old.x <- sort(rep(1:n.nodes, n.groups)) - 1
  old.y <- rep(1:n.groups, n.nodes) - 1
  new.x <- old.x %% width
  new.y <- old.y + floor(old.x/width)*(n.groups+2)
  colblock <- rgb(1-membership.share, 1-membership.share, 1-membership.share)

  plot(range(c(new.x, new.x+1)), -range(c(-1, new.y, new.y+1)), ty="n", axes=FALSE, xlab="", ylab="", main=main)  #, ...)
  rect(new.x, -new.y, new.x+1, -(new.y+1), col=colblock)

  text.x <- (1:n.nodes - 1) %% width + 0.5
  text.y <- floor((1:n.nodes - 1)/width)*(n.groups+2) - 0.5
  text(text.x, -text.y, node.labels)

  #plot the block matrix.
  x.points <- max(new.x) + (-n.groups+1):0
  y.points <- max(new.y) + (-n.groups+1):0

  x.grid <- sort(rep(x.points, n.groups))
  y.grid <- rep(y.points, n.groups)

  #red/white/blue color scheme. Maxes out at a deflection of 6 from the main.
  block.col <- color.block.one(block.matrix)
  rect(x.grid, -y.grid, x.grid+1, -(y.grid+1), col=block.col)
  text(rep(max(new.x - n.groups + 0.5), n.groups), -(y.grid + 0.5), 1:n.groups)
  text(x.points+0.5, -rep(max(new.y - n.groups + 0.5), n.groups), 1:n.groups)

}


single.membership.plot <- function (number.membership,
                                    block.matrix,
                                    main="Stochastic Block Model Outcome",
                                    ...) {
  #number.membership=sample(4, 100, replace=TRUE)

  block.membership <- number.to.vector (number.membership, dim(block.matrix)[1])
  block.membership.plot  (block.membership, block.matrix, main=main, ...)

}


###############################################################
# Coefficient plots.

# modified dot chart, top to bottom.
dotchart.coef <- function (values,
                           names=1:length(values),
                           sd=NULL,
                           interval=NULL,   #n by 2.
                           main="",
                           ...) {
  #values=rnorm(10); names=1:length(values); sd=NULL; interval=NULL

  xlims <- values; if (!is.null(interval)) xlims <- c(interval) else if (!is.null(sd)) xlims <- c(xlims + 2*sd, xlims - 2*sd)

  ypts <- -(1:length(values))
  plot (range(c(xlims, 0)), c(-length(values)-0.5, -0.5), ty="n", axes=FALSE, xlab="Value", ylab="Covariate", main=main, ...)
  box()
  axis(1)
  axis(2, ypts, names, las=2) #, cex.axis=3/sqrt(length(values)))
  abline (h=ypts, col=8, lty=3)
  abline (v=0)

  points (values, ypts, pch=19, col=2, cex=2)
  text (rep(mean(range(c(xlims, 0))), length(values)),
        ypts+0.5,
        names,
        cex=2/sqrt(length(values)))  #added 6-25-14.

  if (!is.null(sd)) {
    if (length(sd) != length(values)) stop ("Length of SD vector does not match length of values.")
    segments(values-2*sd, ypts, values+2*sd, col=2, lwd=3)
  }
  if (!is.null(interval)) {
    if (nrow(interval) != length(values)) stop ("Row count of interval matrix does not match length of values.")
    segments(interval[,1], ypts, interval[,2], col=2, lwd=3)
  }

}


dotchart.two <- function (a.values,
                          b.values,

                          node.names=1:length(a.values),

                          sd.a=NULL,
                          interval.a=NULL,

                          sd.b=NULL,
                          interval.b=NULL,

                          xlab="Sender Effect",
                          ylab="Receiver Effect",

                          main="",
                          ...) {

  xlims <- a.values
  if (!is.null(interval.a)) xlims <- c(interval.a) else if (!is.null(sd.a)) xlims <- c(xlims + 2*sd.a, xlims - 2*sd.a)
  ylims <- b.values
  if (!is.null(interval.b)) ylims <- c(interval.b) else if (!is.null(sd.b)) ylims <- c(ylims + 2*sd.b, ylims - 2*sd.b)

  plot (range(c(xlims, 0)), range(c(ylims, 0)), ty="n", xlab=xlab, ylab=ylab, main=main, ...)

  if (!is.null(sd.a)) {
    if (length(sd.a) != length(a.values)) stop ("Length of SD vector does not match length of a.values.")
    segments(a.values-2*sd.a, b.values, a.values+2*sd.a, b.values, col=2, lwd=3)
    if (length(sd.b) != length(b.values)) stop ("Length of SD vector does not match length of b.values.")
    segments(a.values, b.values-2*sd.b, a.values, b.values+2*sd.b, col=2, lwd=3)
  }
  if (!is.null(interval.a)) {
    if (nrow(interval.a) != length(a.values)) stop ("Row count of interval.a matrix does not match length of a.values.")
    segments(interval.a[,1], b.values, interval.a[,2], b.values, col=2, lwd=3)
    if (nrow(interval.b) != length(b.values)) stop ("Row count of interval.b matrix does not match length of b.values.")
    segments(a.values, interval.b[,1], a.values, interval.b[,2], col=2, lwd=3)

  }

  text (a.values, b.values, node.names)

}





###############################################################
# Modified dendrograms with internal nodes.


simple.split <- function (input) {
  if (length(input)>1) {
    l1 <- sample(c(floor(length(input)/2), ceiling(length(input)/2)))
    output <- list(input[1:l1[1]], input[l1[1]+1:l1[2]])
  } else if (rbinom(1,1,0.5)>0) output <- list(c(), input) else output <- list(input, c())
  return(output)
}

dendrogram.internals <- function (leaf.parents, internal.parents) { #, block.values=rep(0, length(internal.parents)) {
  #internal.parents=c(0,1,1,2,2,1,4,4,1); leaf.parents=sample(9, 100, replace=TRUE)

  #start with ordering.
  n.nodes <- length(leaf.parents)
  done <- 0*internal.parents
  depth <- rep(NA, length(internal.parents))
  #descendant.count <- 0*internal.parents

  current.internal <- which(internal.parents == 0)
  layout.vector <- current.internal
  depth[current.internal] <- 0

  while (any(done==0)) {

    descendants <- which(internal.parents == current.internal)
    depth[descendants] <- depth[current.internal]+1

    pos1 <- which(layout.vector == current.internal)
    if (pos1 > 1) lv.left <- layout.vector[1:(pos1-1)] else lv.left <- NULL
    if (pos1 < length(layout.vector)) lv.right <- layout.vector[(pos1+1):length(layout.vector)] else lv.right <- NULL
    ss1 <- simple.split(descendants)
    layout.vector <- c(lv.left, ss1[[1]], current.internal, ss1[[2]], lv.right)
    if (length(unique(layout.vector)) != length(layout.vector)) stop()

    done[current.internal] <- 1
    if (any(done==0)) {
      current.internal <- which(!is.na(depth) & done == 0)
      if (length(current.internal)>1) current.internal <- sample(current.internal, 1)
    }
  }

  leaf.pos <- new.parent <- c()
  leaf.count <- internal.x.pos <- rep(0, length(internal.parents))
  for (kk in layout.vector) {
    old.length <- length(leaf.pos)
    leaf.pos <- c(leaf.pos, which(leaf.parents == kk), 0, 0)
    new.parent <- c(new.parent, rep(kk, sum(leaf.parents == kk)))
    leaf.count[kk] <- sum(leaf.parents == kk)
    internal.x.pos[kk] <- (length(leaf.pos) + old.length - 3)/2
  }
  actual.x <- (match(1:n.nodes, leaf.pos)-1)/length(leaf.pos)
  internal.x.pos <- internal.x.pos*n.nodes/length(leaf.pos)

  descendant.count <- leaf.count
  for (dd in max(depth):1) {
    depth.set <- which(depth==dd)
    for (ee in depth.set) descendant.count[internal.parents[ee]] <- descendant.count[internal.parents[ee]] + descendant.count[ee]
  }

  return(list(leaf.pos=leaf.pos,
              leaf.count=leaf.count,
              actual.x=actual.x,
              depth=depth,
              new.parent=new.parent,
              descendant.count=descendant.count,
              internal.x.pos=internal.x.pos))
}

circular.dendrogram <- function (leaf.parents,
                                 internal.parents,
                                 block.values=rep(0, length(internal.parents)),
                                 node.labels=1:length(leaf.parents),
                                 ...) {
  #internal.parents=c(0,1,1,2,2,1,4,4,1); leaf.parents=sample(9, 100, replace=TRUE); block.values=rep(0, length(internal.parents))
  n.nodes <- length(leaf.parents)
  n.groups <- length(internal.parents)

  dend.ints <- dendrogram.internals(leaf.parents, internal.parents)
  polar.to.xy <- function(radius, angle)    #angle goes from 0:1.
    radius*cbind(cos(angle*2*pi), sin(angle*2*pi))

  xy.internals <- polar.to.xy (dend.ints$depth/max(dend.ints$depth+1),
                               dend.ints$internal.x.pos/n.nodes)
  xy.internals.2 <- polar.to.xy (dend.ints$depth/max(dend.ints$depth+1) - 0.12,
                                 dend.ints$internal.x.pos/n.nodes)
  xy.leaves <- polar.to.xy (rep(1, length(n.nodes)),
                            dend.ints$actual.x)
  plot (xy.leaves, ty="n", xlab="", ylab="", axes=FALSE, ...)

  #grey rings.
  for (dd in 1:max(dend.ints$depth)) {
    xyr <- polar.to.xy (rep(dd/max(dend.ints$depth+1), 1000),
                        seq(0,1.001, length=1000))
    lines(xyr, col=8, lty=3)
  }

  #Go Leaves Go!
  text (xy.leaves[,1], xy.leaves[,2], node.labels, cex=sqrt(50)/sqrt(n.nodes))

  #arcs.
  for (kk in 1:n.groups) if (sum(leaf.parents==kk) > 0) {
    picks <- range(dend.ints$actual.x[leaf.parents == kk])*n.nodes
    range.pts <- seq(min(picks), max(picks), by=0.01)
    lines (polar.to.xy(rep(0.93, length(range.pts)), range.pts/n.nodes), lwd=2)
  }
  #pointies.
  ins <- polar.to.xy (rep(0.93, n.nodes), dend.ints$actual.x)
  outs <- polar.to.xy (rep(0.97, n.nodes), dend.ints$actual.x)
  segments (ins[,1], ins[,2], outs[,1], outs[,2], lwd=2)


  #internal nodes.
  #arc-ties.
  xy.int.ends <- polar.to.xy (rep(0.93, n.groups), dend.ints$internal.x.pos/n.nodes)
  segments (xy.internals[,1], xy.internals[,2],
            xy.int.ends[,1], xy.int.ends[,2],
            lwd=2)

  #internal ties: just straight lines sucks!

  parent.set <- unique(internal.parents[internal.parents != 0])
  for (kk in parent.set) {
    picks <- dend.ints$internal.x.pos[c(which(internal.parents == kk), kk)]
    subradius <- mean(dend.ints$depth[internal.parents == kk])/(max(dend.ints$depth)+1) - 0.12
    range.pts <- seq(min(picks), max(picks), by=0.01)
    lines (polar.to.xy(rep(subradius, length(range.pts)), range.pts/n.nodes), lwd=2)
  }

  #short lines from master arcs.
  segments (xy.internals[internal.parents != 0, 1], xy.internals[internal.parents != 0, 2],
            xy.internals.2[internal.parents != 0, 1], xy.internals.2[internal.parents != 0, 2],
            lwd=2)



  points(xy.internals,
         pch=19,
         cex=6.2*sqrt(5)/nrow(xy.internals))
  points(xy.internals,
         pch=19,
         cex=6*sqrt(5)/nrow(xy.internals),
         col=color.block.one(block.values))


}



#circular.dendrogram (leaf.parents=sample(9, 100, replace=TRUE), internal.parents=c(0,1,1,2,2,1,4,4,1))
