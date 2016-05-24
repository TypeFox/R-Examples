## This is a patched version of drop.tip that will keep the nodes in
## the correct order.  Otherwise it is exactly the same as the
## drop.tip in ape version 2.4-1 (right down to the use of dim(x)[1]
## instead of nrow(x)).  See REMOVE/REPLACE/DONE at the end for changes.
drop.tip.fixed <- function(phy, tip, trim.internal = TRUE, subtree =
                           FALSE, root.edge = 0, rooted = is.rooted(phy)) {
  if (!inherits(phy, "phylo")) 
    stop("object \"phy\" is not of class \"phylo\"")
  Ntip <- length(phy$tip.label)
  if (is.character(tip)) 
    tip <- which(phy$tip.label %in% tip)
  if (!rooted && subtree) {
    phy <- root(phy, (1:Ntip)[-tip][1])
    root.edge <- 0
  }
  phy <- reorder(phy)
  NEWROOT <- ROOT <- Ntip + 1
  Nnode <- phy$Nnode
  Nedge <- dim(phy$edge)[1]
  if (subtree) {
    trim.internal <- TRUE
    tr <- reorder(phy, "pruningwise")
    N <- node.depth(phy)
  }
  wbl <- !is.null(phy$edge.length)
  edge1 <- phy$edge[, 1]
  edge2 <- phy$edge[, 2]
  keep <- !logical(Nedge)
  if (is.character(tip)) 
    tip <- which(phy$tip.label %in% tip)
  if (!rooted && subtree) {
    phy <- root(phy, (1:Ntip)[-tip][1])
    root.edge <- 0
  }
  keep[match(tip, edge2)] <- FALSE
  if (trim.internal) {
    ints <- edge2 > Ntip
    repeat {
      sel <- !(edge2 %in% edge1[keep]) & ints & keep
      if (!sum(sel)) 
        break
      keep[sel] <- FALSE
    }
    if (subtree) {
      subt <- edge1 %in% edge1[keep] & edge1 %in% edge1[!keep]
      keep[subt] <- TRUE
    }
    if (root.edge && wbl) {
      degree <- tabulate(edge1[keep])
      if (degree[ROOT] == 1) {
        j <- integer(0)
        repeat {
          i <- which(edge1 == NEWROOT & keep)
          j <- c(i, j)
          NEWROOT <- edge2[i]
          degree <- tabulate(edge1[keep])
          if (degree[NEWROOT] > 1) 
            break
        }
        keep[j] <- FALSE
        if (length(j) > root.edge) 
          j <- 1:root.edge
        NewRootEdge <- sum(phy$edge.length[j])
        if (length(j) < root.edge && !is.null(phy$root.edge)) 
          NewRootEdge <- NewRootEdge + phy$root.edge
        phy$root.edge <- NewRootEdge
      }
    }
  }
  if (!root.edge) 
    phy$root.edge <- NULL
  phy$edge <- phy$edge[keep, ]
  if (wbl) 
    phy$edge.length <- phy$edge.length[keep]
  TERMS <- !(phy$edge[, 2] %in% phy$edge[, 1])
  oldNo.ofNewTips <- phy$edge[TERMS, 2]
  n <- length(oldNo.ofNewTips)
  phy$edge[TERMS, 2] <- rank(phy$edge[TERMS, 2])
  if (subtree || !trim.internal) {
    tips.kept <- oldNo.ofNewTips <= Ntip & !(oldNo.ofNewTips %in% 
                   tip)
    new.tip.label <- character(n)
    new.tip.label[tips.kept] <- phy$tip.label[-tip]
    node2tip <- oldNo.ofNewTips[!tips.kept]
    new.tip.label[!tips.kept] <- if (subtree) {
      paste("[", N[node2tip], "_tips]", sep = "")
    }
    else {
      if (is.null(phy$node.label)) 
        rep("NA", length(node2tip))
      else phy$node.label[node2tip - Ntip]
    }
    if (!is.null(phy$node.label)) 
      phy$node.label <- phy$node.label[-(node2tip - Ntip)]
    phy$tip.label <- new.tip.label
  }
  else phy$tip.label <- phy$tip.label[-tip]
  if (!is.null(phy$node.label)) 
    phy$node.label <- phy$node.label[sort(unique(phy$edge[, 
                                                          1])) - Ntip]
  phy$Nnode <- dim(phy$edge)[1] - n + 1L

  ## REMOVE:
  ##     newNb <- integer(n + phy$Nnode)
  ##     newNb[NEWROOT] <- n + 1L
  ##     sndcol <- phy$edge[, 2] > n
  ##     phy$edge[sndcol, 2] <- newNb[phy$edge[sndcol, 2]] <- (n + 
  ##         2):(n + phy$Nnode)
  ##     phy$edge[, 1] <- newNb[phy$edge[, 1]]
  ## REPLACE:
  i <- phy$edge > n
  phy$edge[i] <- match(phy$edge[i], sort(unique(phy$edge[i]))) + n
  ## DONE:
  
  storage.mode(phy$edge) <- "integer"
  collapse.singles(phy)
}

## However, this version is much simpler:
function(phy, tip, trim.internal = TRUE, subtree = FALSE,
         root.edge = 0, rooted = is.rooted(phy)) {
  Ntip <- length(phy$tip.label)
  if (is.character(tip)) 
    tip <- which(phy$tip.label %in% tip)

  phy <- reorder(phy)
  NEWROOT <- ROOT <- Ntip + 1
  Nnode <- phy$Nnode
  Nedge <- nrow(phy$edge)

  wbl <- !is.null(phy$edge.length)
  edge1 <- phy$edge[, 1]
  edge2 <- phy$edge[, 2]
  keep <- !(edge2 %in% tip)  

  ints <- edge2 > Ntip
  repeat {
    sel <- !(edge2 %in% edge1[keep]) & ints & keep
    if (!sum(sel)) 
      break
    keep[sel] <- FALSE
  }

  phy2 <- phy
  phy2$edge <- phy2$edge[keep, ]
  if (wbl) 
    phy2$edge.length <- phy2$edge.length[keep]
  TERMS <- !(phy2$edge[, 2] %in% phy2$edge[, 1])
  oldNo.ofNewTips <- phy2$edge[TERMS, 2]
  n <- length(oldNo.ofNewTips)
  idx.old <- phy2$edge[TERMS, 2]
  phy2$edge[TERMS, 2] <- rank(phy2$edge[TERMS, 2])
  phy2$tip.label <- phy2$tip.label[-tip]
  if (!is.null(phy2$node.label))
    phy2$node.label <-
      phy2$node.label[sort(unique(phy2$edge[, 1])) - Ntip]
  phy2$Nnode <- nrow(phy2$edge) - n + 1L
  i <- phy2$edge > n
  phy2$edge[i] <- match(phy2$edge[i], sort(unique(phy2$edge[i]))) + n
  storage.mode(phy2$edge) <- "integer"
  collapse.singles(phy2)
}
