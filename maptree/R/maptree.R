# maptree package
#   for graphing and mapping of hierarchical clustering and
#   regression trees
# denis white, us epa, 15 November 2006, version 1.4-3
#
# function calls
#
############################################################
#
# clip.clust <- function (cluster, data=NULL, k=NULL, h=NULL)
#
# clip.rpart <- function (tree, cp=NULL, best=NULL) 
#
# draw.clust <- function (cluster, data=NULL, cex=par("cex"), 
#     pch=par("pch"), size=2.5*cex, col=NULL, nodeinfo=FALSE, 
#     membership=FALSE, cases="obs", new=TRUE)
#
# draw.tree <- function (tree, cex=par("cex"), pch=par("pch"), 
#     size=2.5*cex, col=NULL, nodeinfo=FALSE, units="",  
#     cases="obs", digits=getOption("digits"), 
#     print.levels=TRUE, new=TRUE)
#
# group.clust <- function (cluster, k=NULL, h=NULL)
#
# group.tree <- function (tree)
#
# kgs <- function (cluster, diss, alpha=1, maxclust=NULL)
#
# map.groups <- function (pts, group, pch=par("pch"), size=2, 
#     col=NULL, border=NULL, new=TRUE)
#
# map.key <- function (x, y, labels=NULL, cex=par("cex"), 
#    par=par("pch"), size=2.5*cex, col=NULL, head="", 
#    sep=0.25*cex, new=FALSE)
#
# ngon <- function (xydc, n=4, angle=0, type=1)
#
# twins.to.hclust <- function (cluster)


############################################################

clip.clust <- function (cluster, data = NULL, k=NULL, h=NULL)
  # analogous to prune.tree
  # cluster is class hclust or twins
  # data is clustered dataset (provided by twins but not hclust)
  # k is desired number of groups
  # h is height at which to cut for grouping
  # needs k or h, k takes precedence
  # returns pruned cluster
{
  if (is.null (h) && is.null (k))
    stop ("clip.clust: both h=NULL, k=NULL")
  if (!is.null (h) && h > max (cluster$height)) 
    stop ("clip.clust: h > max (height)")
  if (!is.null (k) && (k == 1 || k > length (cluster$height)))
    stop("clip.clust: k==1 || k=>nobs")
  if ("hclust" %in% class (cluster)) 
  {
    if (is.null (data))
      stop ("clip.clust: no data provided for hclust object")
    clust <- cluster
  }
  else if (inherits (cluster, "twins"))
  {
    if (! ("data" %in% names (cluster)))
      if (is.null (data)) 
        stop ("clip.clust: no data provided for twins object")
      else cluster$data <- data
    clust <- twins.to.hclust (cluster)
    # clust <- as.hclust (cluster)
  }
  else
    stop("clip.clust: input not hclust or twins")
  merg <- clust$merge
  hite <- clust$height
  nmerg <- nrow(merg)
  if (!is.null (k)) keep <- rev (order (hite))[1:(k-1)]
  else keep <- seq (nmerg)[hite > h]
  numerg <- matrix (0, nrow=length (keep), ncol=2)
  nuhite <- rep (0, length (keep))

  leaf <- 0
  node <- 0
  trim.clust <- function (oldnode)
  {
    a <- merg[oldnode,1]
    b <- merg[oldnode,2]
    if (match (a, keep, 0) != 0) l <- trim.clust (a)
    else {
      leaf <<- leaf + 1
      l <- -leaf
      }
    if (match (b, keep, 0) != 0) r <- trim.clust (b)
    else {
      leaf <<- leaf + 1
      r <- -leaf
      }
    node <<- node + 1
    numerg[node,] <<- c(l,r)
    nuhite[node] <<- hite[oldnode]
    return (node)
  }

  trim.clust (length (hite))
#  trim.clust(match(max(hite),hite))

  numerg <- matrix (as.integer(numerg), nrow=length (numerg)/2, ncol=2)
  nuhite <- as.double (nuhite - min (nuhite))
  nuordr <- as.double (seq (nrow (numerg) + 1))
  nulabl <- as.character (nuordr)
  g <- group.clust (clust, k, h)
  if ("twins" %in% class(cluster)) data <- clust$data
  m <- split (rownames (data), as.factor (g))
  l <- list (merge=numerg, height=nuhite, order=nuordr, labels=nulabl,
    method=clust$method, call=clust$call, dist.method=clust$dist.method,
    size=table (g), membership=m, data=data)
  class(l) <- class(clust)
  l
}

############################################################

clip.rpart <- function (tree, cp=NULL, best=NULL) 
  # modification to original prune.rpart to add best
{
  ff <- tree$frame
  id <- as.integer (row.names (ff))
  if (is.null (cp)) {
    m <- tree$cptable[, "nsplit"]
    m <- max (m[m < best])
    m <- match (m, tree$cptable[, "nsplit"])
    cp <- tree$cptable[m, "CP"]
    }
  toss <- id[ff$complexity <= cp & ff$var != "<leaf>"]
  if (length (toss) == 0) return(tree)
  newx <- snip.rpart (tree, toss)
  temp <- pmax(tree$cptable[, 1], cp)
  keep <- match (unique (temp), temp)
  newx$cptable <- tree$cptable[keep, ]
  newx$cptable[max(keep), 1] <- cp
  newx
}

############################################################

draw.clust <- function (cluster, data=NULL, cex=par("cex"), 
  pch=par("pch"), size=2.5*cex, col=NULL, nodeinfo=FALSE, 
  cases="obs", new=TRUE)
  # cluster is class hclust or twins
  # data is clustered dataset (provided by twins but not hclust)
  # cex is par parameter, size of text
  # pch is par parameter, shape of symbol at leaves of tree
  # size is in cex units for symbol at leaves of tree
  # if col is NULL, use rainbow()
  # if nodeinfo==TRUE, add a line at each leaf with number
  #   of observations included in leaf
  # cases are names for cluster objects
  # if new=TRUE, call plot.new()
  # returned value is col or generated colors
{
  if ("hclust" %in% class (cluster)) 
  {
    if (is.null (data)) 
      if ("data" %in% names (cluster))
        data <- cluster$data
      else
        stop ("draw.clust: no data provided for hclust object")
    clust <- cluster
  }
  else if (inherits (cluster, "twins")) 
  {
    if (! ("data" %in% names (cluster)))
      if (is.null (data)) 
        stop ("draw.clust: no data provided for twins object")
      else cluster$data <- data
    clust <- twins.to.hclust (cluster)
    # clust <- as.hclust (cluster)
    data <- clust$data
  }
  else stop("draw.clust: input not hclust or twins")
  merg <- clust$merge
  nmerg <- nrow (merg)
  if (nmerg<2) stop ("draw: < 3 clusters")
  if (new) plot.new ()
  hite <- clust$height
  cord <- order (clust$order)
  xmax <- nrow (merg) + 1
  ymax <- max (hite)
  pinx <- par ("pin")[1]
  piny <- par ("pin")[2]
  xmin <- 1
  box <- size * par("cin")[1]
  xscale <- (xmax-xmin)/pinx
  xbh <- xscale * box / 2
  tail <- 0.2
  yscale <- ymax/(piny - tail)
  ytail <- yscale * tail
  ybx <- yscale * box
  ymin <- - ytail
  xf <- 0.1 * (xmax-xmin)
  yf <- 0.1 * (ymax-ymin)
  x1 <- xmin - xf
  x2 <- xmax + xf
  y1 <- ymin - yf
  y2 <- ymax + yf
  par (usr=c(x1, x2, (y1-nodeinfo*ybx), y2))
  if (is.null(col)) kol <- rainbow (xmax)
  else kol <- col
  xmean <- rep (0, nmerg)
  i <- 1
  while (any(xmean == 0)) {
    if (xmean[i] == 0) {
      a <- merg[i,1]
      b <- merg[i,2]
      if (a<0) x1 <- cord[-a] else x1 <- xmean[a]
      if (b<0) x2 <- cord[-b] else x2 <- xmean[b]
      if (x1 != 0 && x2 != 0) xmean[i] <- mean(c(x1,x2))
      }
    i <- i + 1
    if (i > nmerg) i <- 1
    }
  for (i in 1:nmerg) {
    a <- merg[i,1]
    b <- merg[i,2]
    y2 <- hite[i]
    if (a > 0) {
      x1 <- xmean[a]
      y1a <- hite[a]
      }
    else {
      x1 <- cord[-a]
      y1a <- y2 - ytail
      points (x1, y1a-ybx/2, pch=pch, cex=size, col=kol[x1])
      text.default (x1, y1a-(ybx/2), as.character(-a), cex=cex)
      if (nodeinfo) {
        string <- paste (as.character (clust$size[-a]), cases)
        text.default (x1, y1a-1.3*ybx, string, cex=cex)
        }
      }
    if (b > 0) {
      x2 <- xmean[b]
      y1b <- hite[b]
      }
    else {
      x2 <- cord[-b]
      y1b <- y2 - ytail
      points (x2, y1b-ybx/2, pch=pch, cex=size, col=kol[x2])
      text.default (x2, y1b-(ybx/2), as.character(-b), cex=cex)
      if (nodeinfo) {
        string <- paste (as.character (clust$size[-b]), cases)
        text.default (x2, y1b-1.3*ybx, string, cex=cex)
        }
      }
    lines (c(x1,x2),c(y2,y2))
    lines (c(x1,x1),c(y1a,y2))
    lines (c(x2,x2),c(y1b,y2))
    }
  invisible (kol)
}

############################################################

draw.tree <- function (tree, cex=par("cex"), pch=par("pch"),
  size=2.5*cex, col=NULL, nodeinfo=FALSE, units="", 
  cases="obs", digits=getOption("digits"), print.levels=TRUE,
  new=TRUE)
  # tree is object of class tree
  # cex is par parameter, size of text
  # pch is par parameter, shape of symbol at leaves of tree
  # if size=0, draw terminal symbol at leaves 
  #   else symbol with size in cex units
  # if col is NULL, use rainbow()
  # if nodeinfo=TRUE, add a line at each node with mean value
  #   of response, number of observations, and percent
  #   deviance explained (or classified correct)
  # units are for mean value of response (if regression tree)
  # cases are names for observations
  # digits are rounding for response
  # if print.levels=FALSE, do not show factor levels at splits
  # if new=TRUE, call plot.new()
  # returned value is col or generated colors
{
  if (new) plot.new ()
  rtree <- length (attr (tree, "ylevels")) == 0
  tframe <- tree$frame
  rptree <- length (tframe$complexity) > 0
  node <- as.numeric(row.names(tframe))
  depth <- floor (log (node, base=2) + 1e-07)
  depth <- as.vector (depth - min (depth))
  maxdepth <- max(depth)
  x <-  - depth
  y <- x
  leaves <- tframe$var == "<leaf>"
  x[leaves] <- seq(sum(leaves))
  depth <- split(seq(node)[!leaves], depth[!leaves])
  parent <- match(node %/% 2, node)
  left.child <- match(node * 2, node)
  right.child <- match(node * 2 + 1, node)
  for(i in rev(depth))
    x[i] <- 0.5 * (x[left.child[i]] + x[right.child[i]])
  nleaves <- sum(leaves)
  nnodes <- length(node)
  nodeindex <- which (tframe$var != "<leaf>")
  if (rtree) {
    dev <- tframe$dev
    pcor <- rep (0, nnodes)
    for (i in 1:nnodes)
      if (! leaves[i]) {
        l <- dev[node == (node[i]*2)]
        r <- dev[node == (node[i]*2+1)]
        pcor[i] <- dev[i] - l - r
      }
    pcor <- round (pcor/dev[1],3)*100
    }
  else {
    crate <- rep (0, nnodes)
    trate <- 0
    if (! rptree) {
      for (i in 1:nnodes) {
        yval <- tframe$yval[i]
        string <- paste('tframe$yprob[,"',
          as.character(yval), '"]', sep="")
        crate[i] <- eval(parse(text=string))[i]
        if (leaves[i]) trate <- trate + tframe$n[i] * crate[i]
        }
      }
    else {
      for (i in 1:nnodes) {
        yval <- tframe$yval[i]
        nlv <- floor (ncol (tframe$yval2) / 2)
        index <- rev (order (tframe$yval2[i, 2:(nlv+1)]))[1]
        crate[i] <- tframe$yval2[i, (nlv + 1 + index)]
        if (leaves[i]) trate <- trate + tframe$n[i] * crate[i]
        }
      }
    crate <- round (crate,3)*100
    trate <- round (trate/tframe$n[1],3)*100
    }
  if (is.null(col)) kol <- rainbow (nleaves)
  else kol <- col
  xmax <- max(x)
  xmin <- min(x)
  ymax <- max(y)
  ymin <- min(y)
  pinx <- par ("pin")[1]
  piny <- par ("pin")[2]
  xscale <- (xmax - xmin)/pinx
  box <- size * par("cin")[1]
  if (box == 0) xbh <- xscale * 0.2
  else xbh <- xscale * box/2
  chr <- cex * par("cin")[2]
  tail <- box + chr
  yscale <- (ymax - ymin)/(piny - tail)
  ytail <- yscale * tail
  if (box == 0) ybx <- yscale * 0.2
  else ybx <- yscale * box
  ychr <- yscale * chr
  ymin <- ymin - ytail
  xf <- 0.1 * (xmax - xmin)
  yf <- 0.1 * (ymax - ymin)
  x1 <- xmin - xf
  x2 <- xmax + xf
  y1 <- ymin - yf
  y2 <- ymax + yf
  par (usr=c(x1,x2,y1,y2))
  v <- as.character (tframe$var[1])
  if (rptree) {
    sp <- tree$splits[1, ]
    val <- sp["index"]
    if (sp["ncat"] > 1) {
      r <- sp["index"]
      string <- "attributes(tree)$xlevels$"
      string <- paste (string, v, sep="")
      xl <- eval (parse (text=string))
      lf <- rf <- ""
      for (k in 1:sp["ncat"])
        if (tree$csplit[r, k] == 1) 
          lf <- paste (lf, xl[k], sep=",")
        else
          rf <- paste (rf, xl[k], sep=",")
      if (! print.levels) string <- v
      else string <- paste (lf, "=", v, "=", rf)
      }
    else {
      if (sp["ncat"] < 0) op <- "<>" else op <- "><"
      string <- paste (v, op, val)
      }
    }
  else {
    val <- substring(as.character(tframe$splits[1, 1]), 2)
    string <- paste (as.character(v), "<>", val)
    }
  text.default (x[1], y[1], string, cex=cex)
  if (nodeinfo) {
    n <- tframe$n[1]
    if (rtree) {
      z <- round(tframe$yval[1], digits)
      r <- pcor[1]
      string <- paste (z," ",units,"; ",n," ",cases,"; ",
        r,"%",sep="")
      }
    else {
      z <- attr (tree, "ylevels")[tframe$yval[1]]
      r <- crate[1]
      string <- paste (z,"; ",n," ",cases,"; ",r,"%",
        sep="")
      }
    text.default (x[1], y[1]-ychr, string, cex=cex)
    }
  for (i in 2:nnodes) {
    ytop <- ychr * (as.integer(nodeinfo)+1)
    if (y[i] < y[i-1]) {
      lines(c(x[i-1], x[i]), c(y[i-1]-ytop, y[i-1]-ytop))
      lines(c(x[i], x[i]), c(y[i-1]-ytop, y[i]+ychr))
      }
    else {
      lines(c(x[parent[i]], x[i]), c(y[parent[i]]-ytop, 
        y[parent[i]]-ytop))
      lines(c(x[i], x[i]), c(y[parent[i]]-ytop, 
        y[i]+ychr))
      }
    if(! leaves[i]) {
      v <- as.character (tframe$var[i])
      if (rptree) {
        if (length (tree$ordered) > 1) {
          k <- 1
          for (j in 1:(i-1)) {
            m <- tframe$ncompete[j]
            if (m > 0) k <- k + m + 1
            m <- tframe$nsurrogate[j]
            if (m > 0) k <- k + m
            }
          }
        else k <- match (i, nodeindex[-1]) + 1
        sp <- tree$splits[k, ]
        val <- sp["index"]
        if (sp["ncat"] > 1) {
          r <- sp["index"]
          string <- "attributes(tree)$xlevels$"
          string <- paste (string, v, sep="")
          xl <- eval (parse (text=string))
          lf <- rf <- ""
          for (k in 1:sp["ncat"])
            if (tree$csplit[r, k] == 1) 
              lf <- paste (lf, xl[k], sep=",")
            else
              rf <- paste (rf, xl[k], sep=",")
          if (! print.levels) string <- v
          else string <- paste (lf, "=", v, "=", rf)
          }
        else {
          if (sp["ncat"] < 0) op <- "<>" else op <- "><"
          string <- paste (v, op, val)
          }
        }
      else {
        val <- substring(as.character(tframe$splits[i, 1]), 2)
        string <- paste (as.character(v), "<>", val)
        }
      text.default (x[i], y[i], string, cex=cex)
      if (nodeinfo) {
        n <- tframe$n[i]
        if (rtree) {
          z <- round(tframe$yval[i], digits)
          r <- pcor[i]
          string <- paste (z," ",units,"; ",n," ",cases,"; ",
            r,"%",sep="")
          }
        else {
          z <- attr (tree, "ylevels")[tframe$yval[i]]
          r <- crate[i]
          string <- paste (z,"; ",n," ",cases,"; ",r,"%",
            sep="")
          }
        text.default (x[i], y[i]-ychr, string, cex=cex)
        }
      }
    else {
      if (box == 0) {
        lines (c(x[i], x[i]), c(y[i], y[i]+ychr))
        lines (c(x[i]-xbh, x[i]+xbh), c(y[i], y[i]))
        }
      else {
        points (x[i], y[i], pch=pch, cex=size, col=kol[x[i]])
        }
      if (rtree) {
        z <- round(tframe$yval[i], digits)
        text.default(x[i], y[i]-ybx, paste(z,units,sep=" "), cex=cex)
        }
      else {
        z <- attr (tree, "ylevels")[tframe$yval[i]]
        text.default(x[i], y[i]-ybx, z, cex=cex)
        }
      n <- tframe$n[i]
      text.default(x[i], y[i]-ybx-ychr, paste(n,cases,sep=" "), cex=cex)
      if (box != 0)
        text.default (x[i], y[i], as.character(x[i]), cex=cex)
      }
    }
  if (nodeinfo) {
    if (rtree) string <- paste("Total deviance explained =",
      sum(pcor),"%")
    else string <- paste("Total classified correct =",trate,"%")
    if (box == 0) text.default (mean(x),ymin-3*ychr,string,cex=1.2*cex)
    else text.default (mean(x),ymin-1.2*ybx,string,cex=1.2*cex)
    }
}

############################################################

group.clust <- function (cluster, k=NULL, h=NULL)
  # alternative to cutree that orders groups from left to
  # right in draw order
  # cluster is class hclust or twins
  # k is desired number of groups
  # h is height at which to cut for grouping
  # needs k or h, k takes precedence
  # returns vector of membership
{
  if (is.null (h) && is.null (k)) return (cluster$order)
  if (!is.null (h) && h > max (cluster$height)) 
    stop("group.clust: h > max (height)")
  if (!is.null (k) && (k == 1 || k > length (cluster$height)))
    stop("group.clust: k == 1 || k => nobs")
  if ("hclust" %in% class (cluster)) clust <- cluster
  else if (inherits (cluster, "twins"))
    clust <- as.hclust (cluster)
  else
    stop("group.clust: input not hclust or twins")
  merg <- clust$merge
  hite <- clust$height
  ordr <- clust$order
  nmerg <- nrow (merg)
  group <- rep (0, nmerg+1)
  if (!is.null (k)) keep <- rev (order (hite))[1:(k-1)]
  else keep <- seq (nmerg)[hite > h]

  mark.group <- function (node, grup)
  {
    a <- merg[node,1]
    b <- merg[node,2]
    if (a < 0) group[-a] <<- grup
    else mark.group (a, grup)
    if (b < 0) group[-b] <<- grup
    else mark.group (b, grup)
    invisible()
  }

  grup <- 0
  find.group <- function (node)
  {
    a <- merg[node,1]
    b <- merg[node,2]
    if (match (a, keep, 0) != 0) find.group (a)
    else {
      grup <<- grup + 1
      if (a > 0) mark.group (a, grup)
      else group[-a] <<- grup
      }
    if (match (b, keep, 0) != 0) find.group (b)
    else {
      grup <<- grup + 1
      if (b > 0) mark.group (b, grup)
      else group[-b] <<- grup
      }
    invisible ()
  }

  find.group (length (hite))
  grup <- match (grup, unique (grup[clust$order]))
  invisible (group)
}

############################################################

group.tree <- function (tree)
  # alternative to tree$where that orders groups from left
  # to right in draw order
{
  group <- match (tree$where, sort (unique (tree$where)))
  names (group) <- names (tree$where)
  invisible (group)
}

############################################################

kgs <- function (cluster, diss, alpha=1, maxclust=NULL)
  # cluster is class hclust or twins
  # diss is class dist or dissimilarity
  # alpha is weight for number of clusters
  # maxclust is maximum number of clusters to compute for;
  #   if NULL, use n-1
  # needs {maptree}
  # this implementation of complexity O(n*n*maxclust);
  #   needs memory from level to level to cut down compares
  # ref: Kelley LA, Gardner SP, Sutcliffe MJ. 1996. An
  #   automated approach for clustering an ensemble of
  #   NMR-derived protein structures into conformationally-
  #   reated subfamilies. Protein Engineering 9:1063-1065.
{
    spread <- function (mem, diss)
    {
      if (length (mem) > 2) comb <- combn (mem, 2)
      else comb <- matrix (mem, nrow=2, ncol=1)
      n <- ceiling (sqrt (2 * length (diss)))
      sp <- 0
      for (k in seq (ncol (comb))) {
        i <- comb[1, k]
        j <- comb[2, k]
        sp <- sp + diss[n*(i-1) - i*(i-1)/2 + j-i]
      }
      sp <- sp * 2 / (n * (n-1))
      sp
    }
  if ("hclust" %in% class (cluster)) clust <- cluster
  else if (inherits (cluster, "twins"))
    clust <- as.hclust (cluster)
  else
    stop("kgs: input not hclust or twins")
  if (class (diss) != "dist" && 
      class (diss) != "dissimilarity") 
    stop ("kgs: input not dist or dissimilarity")
  n <- length (clust$order)
  if (is.null (maxclust) || maxclust > (n-1)) m <- n - 1
  else m <- maxclust
  avsp <- rep (0, m-1)
  for (i in 2:m) {
    gl <- group.clust (clust, k=i)
    sz <- table (gl)
    nsp <- sp <- 0
    for (j in seq (i)) {
      if (sz[j] > 1) {
        mem <- seq (n)[gl == j]
        sp <- sp + spread (mem, diss)
        nsp <- nsp + 1
    } }
    avsp[i-1] <- sp / nsp
  }
  avsp <- (m-1)*(avsp - min(avsp))/diff (range (avsp)) + 1
  kgs <- avsp + alpha*(2:m)
  names (kgs) <- as.character (2:m)
  kgs
}

############################################################

map.groups <- function (pts, group, pch=par("pch"), size=2, 
  col=NULL, border=NULL, new=TRUE)
  # pts must have components "x" and "y";
  # group is vector of length of number of cases (either
  #   polygons or points) and indexes colors;
  #   if pts are for polygons, then names(group) must match
  #   with pts$x[is.na(pts$y)]
  # if nrow(pts) != length(group) then map with polygon, 
  #   else if pch < 100 map with points, 
  #   else map with ngon (..., n=pch-100)
  # pch is par parameter, shape of point symbol
  # size is in cex units of point symbol
  # if col is NULL, use rainbow()
  # if border is NULL, use fill colors (col),
  #   else the specified color(s)
  # if new=TRUE, call plot.new()
{
  n <- colnames (pts)
  if (! any (n == "x"))
    stop ("map.groups: pts has no $x")
  else if (! any (n == "y"))
    stop ("map.groups: pts has no $y")
  if (nrow(pts) != length(group)) 
    if (is.null (names (group)))
      stop ("map.groups: group has no names")
  if (new) plot.new ()
  nna <- which (! is.na (pts$y))
  rx <- range (pts$x[nna])
  ry <- range (pts$y[nna])
  dx <- diff (rx)
  dy <- diff (ry)
  ex <- 0.02
  plot.window (rx + c(-ex*dx,ex*dx), 
               ry + c(-ex*dy,ex*dy), asp=1)
  dense <- sort (unique (group))
  nc <- length (dense)
  if (is.null(col)) fkol <- rainbow (nc)
  else fkol <- rep (col, nc)
  if (is.null(border)) bkol <- fkol
  else bkol <- rep (border, nc)
  if (nrow (pts) != length (group)) {
    i <- pts$x[is.na (pts$y)]
    j <- match (i, as.integer (names (group)))
    polygon (pts$x, pts$y, lwd=0.1,
      col=fkol[match (group[j], dense)], 
      border=bkol[match (group[j], dense)])
    }
  else if (pch < 100 | mode(pch) == "character") points (pts$x, 
    pts$y, col=fkol[match (group, dense)], pch=pch, cex=size*1.5)
  else apply (data.frame (x=pts$x, y=pts$y, 
    d=size*25.4*par("cex")*par("cin")[1], 
    c=I(fkol[match (group, dense)])), 
    1, ngon, n=pch-100, type=1)
  invisible (fkol)
}

############################################################

map.key <- function (x, y, labels=NULL, cex=par("cex"), 
  pch=par("pch"), size=2.5*cex, col=NULL, head="", 
  sep=0.25*cex, new=FALSE)
  # x,y are lower left coordinates of key 
  #   in proportional units (0-1)
  # labels is vector of labels for classes, or if NULL,
  #   then integers 1:length(col), or "1"
  # cex is par parameter, size of text
  # if pch < 100, use points for symbol, 
  #   else ngon (..., n=pch-100)
  # size is in cex units for key symbols
  # if col is NULL, use rainbow()
  # cex is par parameter, size of text
  # head is text heading for key
  # sep is separation in cex units between adjacent symbols
  #   if sep=0 assume continuous scale and use gon=4
  #   and put lables at breaks between squares
  # if new=TRUE, call plot
  # returned value is col or generated colors
{
  if (is.null (labels))
    if (is.null (col)) labels <- as.vector ("1")
    else labels <- as.character (seq (length (col) - 1))
  nsym <- length (labels)
  if (sep == 0) nsym <- nsym - 1
  if (is.null (col)) kol <- rainbow (nsym)
  else kol <- col
  if (new)
    plot(c(0,1), c(0,1), type="n", axes=FALSE, xlab="", ylab="")
  oldadj <- par ("adj")
  par (adj=0)
  u <- par ("usr")
  ex <- par ("pin")[1]
  ey <- par ("pin")[2]
  uxr <- u[2] - u[1]
  uyr <- u[4] - u[3]
  halfx <- (size * par("cin")[1]) / 3
  xstep <- halfx + 0.05
  ystep <- (size + sep) * par("cin")[2] / 2.5
  px <- x * uxr + u[1]
  py <- y * uyr + u[3]
  hx <- halfx * uxr / ex
  dx <- xstep * uxr / ex
  dy <- ystep * uyr / ey
  qx <- px
  qy <- py - dy
  if (sep == 0) {
    for (i in 1:nsym) {
      qy <- qy + dy
      points (qx, qy, pch=15, cex=size*1.1, col=kol[i])
      text (qx+dx, qy - dy/2, labels[i], cex=cex) }
    text (qx+dx, qy + dy/2, labels[nsym+1], cex=cex)
    }
  else
    for (i in 1:nsym) {
      qy <- qy + dy
      if (pch < 100  | mode(pch) == "character") points (qx,
        qy, col=kol[i], pch=pch, cex=size)
      else ngon (c(qx, qy, size=15*size*par("cin")[1], col=kol[i]), 
        n=pch-100, type=1)
      text (qx+dx, qy, labels[i], cex=cex) }
  if (length (head) > 0)
    {
    qy <- qy + (dy * length (grep ("$", head)))
    if (sep == 0) qy <- qy + 0.5 * dy
    text (qx-hx, qy, head, cex=cex)
    }
  par (adj=oldadj)
  invisible (kol)
}

############################################################

ngon <- function (xydc, n=4, angle=0, type=1)
  # draw or fill regular polygon
  # xydc a four element vector with
  #   x and y of center, d diameter in mm, and c color
  # n number of sides of polygon, n>8 => circle
  #   if n odd, vertex at (0,d/2), else midpoint of side
  # angle is in degrees by which to rotate the figure
  # type=1 => interior filled, type=2 => edge
  # type=3 => both
{
      # scale factors for d based on n (ignoring angle)
      # n = 3, s = (2 + sqrt(3)) / 4     = 0.9330127
      # n = 4, s = 1 / sqrt(2)           = 0.7071068
      # n = 5, s = (1 + cos(.2*pi)) / 2  = 0.9045085
      # n = 6, s = sqrt(3) / 2           = 0.8660254
  u <- par ("usr")
  p <- par ("pin")
  d <- as.numeric (xydc[3])
  inch <- d / 25.4
  s <- 1
  switch (n, stop ("ngon: n=1"), 
             stop ("ngon: n=2"),
             s <- 0.9330127,
             s <- 0.7071068,
             s <- 0.9045085,
             s <- 0.8660254)
  inch <- inch / s
  rad <- inch*((u[2]-u[1])/p[1])/2
  ys <- inch*((u[4]-u[3])/p[2])/2/rad
  if (n > 8) n <- d*4 + 1
  th <- pi*2/n
  costh <- cos (th)
  sinth <- sin (th)
  x <- y <- rep (0,n+1)
  if (n %% 2) {
    x0 <- 0
    y0 <- rad
    }
  else {
    x0 <- -rad*sin(th/2)
    y0 <-  rad*cos(th/2)
    }
  a <- pi*angle/180
  x[1] <- x0*cos(a) - y0*sin(a)
  y[1] <- x0*sin(a) + y0*cos(a)
  for (i in 2:(n+1)) {
    xl <- x[i-1]
    yl <- y[i-1]
    x[i] <- xl*costh - yl*sinth
    y[i] <- xl*sinth + yl*costh
    }
  x <- x + as.numeric (xydc[1])
  y <- y*ys + as.numeric (xydc[2])
  if (type %% 2) polygon (x, y, col=xydc[4], border=0)
  if (type %/% 2) lines (x, y, col=xydc[4])
  invisible ()
}

############################################################

twins.to.hclust <- function (cluster)
{
  if (! inherits(cluster,"twins"))
    stop ("twins.to.hclust: input not twins")
  merg <- cluster$merge
  hite <- cluster$height
  ordr <- cluster$order
  nuhite <- rep (0, length(hite))

  hite.clust <- function (node)
  {
    a <- merg[node,1]
    b <- merg[node,2]
    if (a < 0) l <- rep(match(-a,ordr),2)
    else l <- hite.clust(a)
    if (b < 0) r <- rep(match(-b,ordr),2)
    else r <- hite.clust(b)
    if (r[1] - l[2] == 1) nuhite[node] <<- hite[l[2]]
    return (c(l[1],r[2]))
  }

  n <- length(hite)
  hite.clust (n)

  l <- list()
  l[[1]] <- merg
  l[[2]] <- nuhite
  l[[3]] <- ordr
  l[[4]] <- cluster$order.lab
  l[[5]] <- "unknown"
  l[[6]] <- attr(cluster,"Call")
  l[[7]] <- "unknown"
  l[[8]] <- cluster$data
  names(l) <- 
    c("merge","height","order","labels","method",
      "call","dist.method","data")
  class(l) <- "hclust"
  l
}
