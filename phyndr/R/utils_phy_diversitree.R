get_descendants <- function(node, tree, tips.only=FALSE,
                            edge.index=FALSE) {
  n_tip <- length(tree$tip.label)
  if (length(node) != 1) {
    stop("'node' must be scalar")
  }
  if (is.character(node)) {
    node <- match(node, tree$node.label) + n_tip
    if (is.na(node)) {
      stop(sprintf("Node '%s' not found in tree"), node)
    }
  } else {
    node <- check_integer(node)
    if (node >= 1 && node < tree$Nnode) { # on 1..(n.node), probably
      node <- node + n_tip
    } else if ( !(node > n_tip && node <= n_tip + tree$Nnode)) {
      stop("Invalid node number")
    }
  }

  edge <- tree$edge

  desc <- descendants_C(node, edge, n_tip)
  if (tips.only) {
    desc <- desc[desc <= n_tip]
  }
  if (edge.index) {
    desc <- match(desc, edge[, 2])
  }
  desc
}

##' @useDynLib phyndr
descendants_C <- function(node, edge, n_tip) {
  storage.mode(edge) <- "integer"
  storage.mode(node) <- "integer"
  storage.mode(n_tip) <- "integer"
  .Call("r_descendants", node, edge, n_tip, PACKAGE="phyndr")
}

get_descendants_names <- function(node, phy) {
  phy$tip.label[get_descendants(node, phy, tips.only=TRUE)]
}

mrca_tipset <- function(phy, tips) {
  if (is.character(tips)) {
    tips <- match(tips, phy$tip.label)
  }

  if (length(tips) == 1L) {
    tips
  } else {
    anc <- ancestors(phy, tips)
    j <- which(apply(anc, 1, function(x) length(unique(x))) > 1L)[1]
    anc[j - 1L, 1]
  }
}

ancestors <- function(phy, i=seq_along(phy$tip.label)) {
  anc <- i
  edge <- phy$edge
  while (any(!is.na(i))) {
    i <- edge[match(i, edge[,2]),1]
    anc <- cbind(anc, i, deparse.level=0)
  }

  apply(anc, 1, function(x) c(rev(x[!is.na(x)]), rep(NA, sum(is.na(x)))))
}

## Check that a number can reasonably be considered an integer.
check_integer <- function(x, name=deparse(substitute(x))) {
  if (is.null(x)) {
    stop("NULL argument for ", name)
  }
  nna <- !is.na(x)
  if (length(x) > 0 && !any(nna)) {
    stop("No non-NA values for ", name)
  }
  if (length(x) && max(abs(x[nna] - round(x[nna]))) > 1e-8) {
    stop("Non-integer argument for ", name)
  }
  storage.mode(x) <- "integer"
  x
}
