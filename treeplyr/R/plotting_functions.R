#' Internal functions for painted_clades
.plotRegimes <- function(tree, col=NULL, lwd=1, pal=rainbow, ...){
  if(is.null(col)){
    regNames <- unique(names(unlist(tree$maps)))
    nreg <- length(regNames)
    col <- setNames(pal(nreg), regNames)
  }
  #nodecols <- col[sapply(tree$maps, function(x) names(x)[1])]
  tmp <- plot(tree, edge.color="#FFFFFF00", use.edge.length=TRUE, ...)
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  #if(lastPP$type != "phylogram") stop("Currently only able to plot phylograms")
  nbranch <- nrow(tree$edge)
  .getBranchCoords <- function(i){
    xx <- lastPP$xx[tree$edge[i,]]
    yy <- lastPP$yy[tree$edge[i,]]
    xdist <- diff(xx)
    ydist <- diff(yy)
    map <- tree$maps[[i]]
    cs <- cumsum(c(0, map))/sum(map)
    colmap <- col[names(map)]
    return(list(xx=xx, yy=yy, xdist=xdist, ydist=ydist, cs=cs, colmap=colmap, nsegs=length(cs)-1, segreg = names(colmap)))
  }
  coords <- lapply(1:nbranch, .getBranchCoords)
  .phylogramLines <- function(x){
    xdist <- x$xdist; ydist <- x$ydist; xx <- x$xx; yy <- x$yy
    cs <- x$cs; nsegs <- x$nsegs; segreg <- x$segreg; colmap <- x$colmap
    if(lastPP$direction %in% c("upwards", "downwards")){
      xcoord <- rbind(xx, matrix(xx[2], nrow=nsegs, ncol=2))
      ycoord <- rbind(rep(yy[1],2), cbind(cs[1:(length(cs)-1)]*ydist+yy[1], cs[2:(length(cs))]*ydist+yy[1]))
      rownames(xcoord) <- rownames(ycoord) <- c(segreg[1], segreg)
      cols <- c(colmap[1], colmap)
      dum <- lapply(1:(nsegs+1), function(i) lines(xcoord[i,], ycoord[i, ], col=cols[i], lwd=lwd))
    }
    if(lastPP$direction %in% c("leftwards", "rightwards")){
      ycoord <- rbind(yy, matrix(yy[2], nrow=nsegs, ncol=2))
      xcoord <- rbind(rep(xx[1],2), cbind(cs[1:(length(cs)-1)]*xdist+xx[1], cs[2:(length(cs))]*xdist+xx[1]))
      rownames(xcoord) <- rownames(ycoord) <- c(segreg[1], segreg)
      cols <- c(colmap[1], colmap)
      dum <- lapply(1:(nsegs+1), function(i) lines(xcoord[i,], ycoord[i, ], col=cols[i], lwd=lwd))
    }
  }
  .cladogramLines <- function(x){
    xdist <- x$xdist; ydist <- x$ydist; xx <- x$xx; yy <- x$yy
    cs <- x$cs; nsegs <- x$nsegs; segreg <- x$segreg; colmap <- x$colmap
    xcoord <- cbind(cs[1:(length(cs)-1)]*xdist+xx[1], cs[2:(length(cs))]*xdist+xx[1])
    ycoord <- cbind(cs[1:(length(cs)-1)]*ydist+yy[1], cs[2:(length(cs))]*ydist+yy[1])
    rownames(xcoord) <- rownames(ycoord) <- segreg
    cols <- colmap
    dum <- lapply(1:nsegs, function(i) lines(xcoord[i,], ycoord[i, ], col=cols[i], lwd=lwd))
  }
  .fanLines <- function(x){
    xdist <- x$xdist; ydist <- x$ydist; xx <- x$xx; yy <- x$yy
    cs <- x$cs; nsegs <- x$nsegs; segreg <- x$segreg; colmap <- x$colmap
    circular.plot(lastPP$edge, lastPP$Ntip, lastPP$Nnode, lastPP$xx, lastPP$yy, )
  }
  if(lastPP$type=="fan") warning("type='fan' not currently supported, plotting a radial cladogram")
  plotfn <- switch(lastPP$type, phylogram=.phylogramLines, cladogram=.cladogramLines, unrooted=.cladogramLines, radial=.cladogramLines, fan=.cladogramLines)
  dum <- lapply(coords, plotfn)
  return(col)
}

#' Internal function taken from bayou
.identifyBranches <- function(tree, n, fixed.loc=TRUE, plot.simmap=TRUE){
  mar.old <- par('mar')
  par(mfrow=c(1,1), mar=c(0.1,0.1,0.1,0.1))
  tree <- reorder(tree,"postorder")
  plot(tree, cex=0.5)
  L <- get("last_plot.phylo",envir=.PlotPhyloEnv)
  nH <- .nodeHeights(tree)
  xx <- L$xx
  yy <- L$yy
  yH <- yy[tree$edge[,2]]
  sb <- numeric()
  loc.x <- numeric()
  for(i in 1:n){
    lx <- locator(1)    
    xmatch <- which(nH[,1] < lx$x & nH[,2] > lx$x)
    ymatch <- which(abs(yH[xmatch]-lx$y)==min(abs(yH[xmatch]-lx$y)))
    sb[i] <- xmatch[ymatch]
    loc.x[i] <- lx$x
    points(lx$x,yH[xmatch][ymatch],pch=21,bg="gray50")
    text(lx$x,yH[xmatch][ymatch],labels=sb[i],pos=4)
  }
  loc <- loc.x - nH[sb,1]
  if(plot.simmap){
    cache <- .prepare.branches(tree)
    pars <- list(sb=sb, loc=loc, t2=1:n+1)
    tr <- .toSimmap(.pars2map(pars, cache),cache)
    cols <- .plotRegimes(tr, lwd=1, pal=rainbow)
    L <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    legend(min(L$xx), max(L$yy), legend=names(cols), lwd=3, col=cols)
  }
  par(mar=mar.old)
  out <- list(sb=sb)
  if(fixed.loc) out$loc <- loc
  return(out)
}

#' Internal function to determine tip regimes:
.tipregime <- function(pars, tree){
  ntips <- length(tree$tip.label)
  tree <- reorder(tree, "postorder")
  cache <- .prepare.branches(tree)
  tr <- .toSimmap(.pars2map(pars, cache),cache)
  tip.reg <- as.numeric(sapply(tr$map, function(x) names(x)[length(x)]))
  tip.reg <- tip.reg[which(tr$edge[,2] <= ntips)]
  o <- order(tr$edge[,2][which(tr$edge[,2] <= ntips)])
  tip.reg <- tip.reg[o]
  tip.reg <- setNames(tip.reg, tr$tip.label)
  return(tip.reg)
}

.toSimmap <- function(map, cache){
  maps <- lapply(1:length(cache$edge.length), function(x){ y <- map$segs[which(map$branch==x)]; names(y) <- map$theta[which(map$branch==x)]; y })  
  tree <- cache$phy
  tree$maps <- maps
  return(tree)
}

#' This is an internal function modified from geiger's function .prepare.bm.univariate for use with OU models.
.prepare.branches <- function(tree, ...){
  ntips <- length(tree$tip.label)
  tree <- reorder(tree, "postorder")
  ## function from geiger 2.0, written by Jonathan Eastman
  .cache.descendants <- function (phy) {
    N = as.integer(Ntip(phy))
    n = as.integer(Nnode(phy))
    phy = reorder(phy, "postorder")
    zz = list(N = N, MAXNODE = N + n, ANC = as.integer(phy$edge[,1]), DES = as.integer(phy$edge[, 2]))
    res = .Call("geiger_descendants", phy = zz, package = "treeplyr")
    return(res)
  }
  desc <- .cache.descendants(tree)
  plook <- function(x){mapply(paste,x[2:length(x)],x[1:(length(x)-1)],sep=",")}
  tB <- desc$anc[1:ntips]
  tB <- mapply(c,1:ntips,tB, SIMPLIFY=FALSE)
  lookup <- lapply(tB,plook)
  edge.names <- mapply(paste, tree$edge[,1], tree$edge[,2],sep=",")
  bdesc <- lapply(edge.names,function(branch) which(edge.names %in% unique(unlist(sapply(lookup,function(look) if(branch %in% look) look[1:which(branch==look)])))))
  bdesc <- lapply(bdesc,function(x) x[-length(x)])
  cache <- list(n=ntips, edge=tree$edge, edge.length=tree$edge.length, bdesc=bdesc, desc=desc, phy=tree)
  #ntips <- length(tree$tip.label)
  #rownames(tree$edge) <- 1:(length(tree$edge[,1]))
  #cache <- .prepare.bm.univariate(tree, X, SE=SE)#, ...)
  #ind <- as.numeric(rownames(cache$edge))
  #cache$n <- ntips
  #cache$N <- nrow(cache$phy$edge)
  #cache$nH <- phytools::nodeHeights(tree)[ind,1]
  #if(is.null(pred)){
  #  pred <- cbind(rep(0, ntips))
  #  rownames(pred) <-  cache$tip.label
  ##}
  #cache$maps <- tree$maps[ind]
  #cache$mapped.edge <- tree$mapped.edge[ind,]
  ##cache$height <- max(phytools::nodeHeights(tree))
  #cache$anc <- cache$phy$edge[,1]
  #cache$des <- cache$phy$edge[,2]
  #cache$ntips <- length(X)
  #cache$ind <- ind
  #cache$ordering <- "postorder"
  ##cache$ht <- .heights.cache(cache)
  #cache$edge <- unname(cache$edge)
  #o <- match(tree$tip.label, cache$tip.label)
  #cache$pred <- cbind(pred[o,])
  #  rownames(cache$edge)=NULL
  ##cache$distFromRoot <- .pruningwise.distFromRoot(cache$phy, n=cache$n, N=cache$N)
  ##cache$tipFromRoot <- cache$distFromRoot[1:cache$n]
  #cache$ultrametric <- as.numeric(is.ultrametric(cache$phy))
  #D <- max(cache$distFromRoot[1:cache$n]) - cache$distFromRoot[1:cache$n]
  #cache$D <- D - mean(D)
  #phy2 <- cache$phy
  #for (i in 1:cache$n) {
  #  tmp = phy2$edge.length[which(cache$phy$edge[,2]==i)]
  #  phy2$edge.length[which(cache$phy$edge[,2]==i)] = tmp + cache$D[i]      
  #}
  #cache$times <- .pruningwise.branching.times(phy2, n=cache$n, des=phy2$edge[,2], anc=phy2$edge[,1])
  #names(cache$times) <- (cache$n+1):(cache$n+cache$phy$Nnode) 
  #cache$Tmax <- max(cache$times)
  #cache$branches.anc <- lapply(1:cache$N, function(x) which(names(cache$times)==cache$phy$edge[,1][x]))
  #cache$branches.des <- lapply(1:cache$N, function(x) which(names(cache$times)==cache$phy$edge[,2][x]))
  #cache$branches.des2 <- cache$branches.des
  #cache$branches.des2[which(sapply(cache$branches.des,length)==0)] <- -1
  #cache$branches.des2 <- unlist(cache$branches.des2)
  #cache$branches.anc2 <- unlist(cache$branches.anc)
  #cache$externalEdge <- cache$phy$edge[,2]<=cache$n
  return(cache)
}

## New version of .pars2map is faster, returns 3 elements rather than 2 named elements
.pars2map <- function(pars, cache){
  nbranch <- length(cache$edge.length)
  nshifts <- table(pars$sb)
  shifts <- rep(0,nbranch)
  shifts[as.numeric(attributes(nshifts)$dimnames[[1]])]<- nshifts
  irow <- rep(1:nbranch,shifts+1)
  segs <- c(cache$edge.length, pars$loc)
  tmp.o <- c(1:nbranch, pars$sb)
  #names(segs) <- tmp.o
  add.o <- order(tmp.o,segs)
  segs <- segs[add.o]
  ind <- tmp.o[add.o]
  #ind <- tmp.o
  t2index <- add.o[which(add.o > nbranch)]
  t2b <- c(rep(1,length(segs)))
  t2b[match(t2index,add.o)+1] <- pars$t2[t2index-nbranch]
  loc.o <- order(pars$loc,decreasing=TRUE)
  sandwiches <- loc.o[which(duplicated(pars$sb[loc.o]))]
  if(length(sandwiches)>0){
    sb.down <- pars$sb[-sandwiches]
    t2.down <- pars$t2[-sandwiches]
  } else {sb.down <- pars$sb; t2.down <- pars$t2}
  sb.o <- order(sb.down)
  sb.down <- sb.down[sb.o]
  t2.down <- t2.down[sb.o]
  sb.desc <- cache$bdesc[sb.down]
  desc.length <- unlist(lapply(sb.desc, length),F,F)
  sb.desc <- sb.desc[which(desc.length>0)]
  #names(t2b) <- names(segs)
  sb.desc2 <- unlist(sb.desc,F,F)
  sb.dup <- duplicated(sb.desc2)
  sb.desc3 <- sb.desc2[which(!sb.dup)]
  t2.names <- rep(t2.down[which(desc.length>0)], unlist(lapply(sb.desc,length),F,F))
  t2.names <- t2.names[which(!sb.dup)]
  #t2b[as.character(unlist(sb.desc3,F,F))] <- t2.names
  t2b[match(sb.desc3, ind)] <- t2.names
  base <- duplicated(ind)*c(0,segs[1:(length(segs)-1)])
  segs <- segs-base
  #maps <- lapply(1:nbranch, function(x) segs[ind==x])
  #maps <- lapply(maps, function(x) if(length(x) >1) {c(x[1],diff(x[1:length(x)]))} else x)
  return(list(segs=segs,theta=t2b, branch=ind))
}

## nodeHeights function from phytools
.nodeHeights <- function (tree) {
  if (!inherits(tree, "phylo")) 
    stop("tree should be an object of class \"phylo\".")
  if (attr(tree, "order") != "cladewise" || is.null(attr(tree, 
                                                         "order"))) 
    t <- reorder(tree)
  else t <- tree
  root <- length(t$tip.label) + 1
  X <- matrix(NA, nrow(t$edge), 2)
  for (i in 1:nrow(t$edge)) {
    if (t$edge[i, 1] == root) {
      X[i, 1] <- 0
      X[i, 2] <- t$edge.length[i]
    }
    else {
      X[i, 1] <- X[match(t$edge[i, 1], t$edge[, 2]), 2]
      X[i, 2] <- X[i, 1] + t$edge.length[i]
    }
  }
  if (attr(tree, "order") != "cladewise" || is.null(attr(tree, 
                                                         "order"))) 
    o <- apply(matrix(tree$edge[, 2]), 1, function(x, y) which(x == 
                                                                 y), y = t$edge[, 2])
  else o <- 1:nrow(t$edge)
  return(X[o, ])
}

## Functions taken from dplyr
group_by_prepare <- function (.data, ..., .dots, add = FALSE) {
  new_groups <- lazyeval::all_dots(.dots, ...)
  is_name <- vapply(new_groups, function(x) is.name(x$expr), 
                    logical(1))
  has_name <- names2(new_groups) != ""
  needs_mutate <- has_name | !is_name
  if (any(needs_mutate)) {
    .data <- mutate_(.data, .dots = new_groups[needs_mutate])
  }
  new_groups <- lazyeval::auto_name(new_groups)
  groups <- lapply(names(new_groups), as.name)
  if (add) {
    groups <- c(groups(.data), groups)
  }
  groups <- groups[!duplicated(groups)]
  list(data = .data, groups = groups)
}

names2 <- function (x) {
  names(x) %||% rep("", length(x))
}

'%||%'<- function (x, y) if (is.null(x)) y else x
