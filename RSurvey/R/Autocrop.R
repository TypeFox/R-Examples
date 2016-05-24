# Approximate the shape of an area defined by a set of points.

Autocrop <- function(mesh, max.len, max.itr=10000) {

  ## Additional functions

  # Return outer elements with arc lengths greater than max.len

  ModTri <- function(tri) {
    n <- nrow(tri)

    arc.id <- rep(NA, n)
    arc.id[tri[, "tr1"] == 0] <- tri[tri[, "tr1"] == 0, "arc1"]
    arc.id[tri[, "tr2"] == 0] <- tri[tri[, "tr2"] == 0, "arc2"]
    arc.id[tri[, "tr3"] == 0] <- tri[tri[, "tr3"] == 0, "arc3"]

    elem.arc <- rep(NA, n)
    elem.arc[arc.id == tri[, "arc1"]] <- 1
    elem.arc[arc.id == tri[, "arc2"]] <- 2
    elem.arc[arc.id == tri[, "arc3"]] <- 3

    pt1 <- array(elem.build[elem.arc, 1])
    pt2 <- array(elem.build[elem.arc, 2])

    pt1.id <- pt2.id <- rep(NA, n)
    pt1.id[pt1 == 1 & !is.na(pt1)] <- tri[pt1 == 1 & !is.na(pt1), "node1"]
    pt1.id[pt1 == 2 & !is.na(pt1)] <- tri[pt1 == 2 & !is.na(pt1), "node2"]
    pt1.id[pt1 == 3 & !is.na(pt1)] <- tri[pt1 == 3 & !is.na(pt1), "node3"]
    pt2.id[pt2 == 1 & !is.na(pt2)] <- tri[pt2 == 1 & !is.na(pt2), "node1"]
    pt2.id[pt2 == 2 & !is.na(pt2)] <- tri[pt2 == 2 & !is.na(pt2), "node2"]
    pt2.id[pt2 == 3 & !is.na(pt2)] <- tri[pt2 == 3 & !is.na(pt2), "node3"]

    pt1xy <- pt2xy <- list(x=rep(NA, n), y=rep(NA, n))
    pt1xy$x[!is.na(pt1.id)] <- mesh$x[na.omit(pt1.id)]
    pt1xy$y[!is.na(pt1.id)] <- mesh$y[na.omit(pt1.id)]
    pt2xy$x[!is.na(pt2.id)] <- mesh$x[na.omit(pt2.id)]
    pt2xy$y[!is.na(pt2.id)] <- mesh$y[na.omit(pt2.id)]

    arcLength <- sqrt((pt2xy$x - pt1xy$x)^2 + (pt2xy$y - pt1xy$y)^2)

    omit.elems <- tri[!is.na(arcLength) & arcLength >  max.len, "elem"]

    new.tri <- NA
    if (length(omit.elems) > 0) {
      tri[tri[, "tr1"] %in% omit.elems, c("tr1", "arc1")] <- 0
      tri[tri[, "tr2"] %in% omit.elems, c("tr2", "arc2")] <- 0
      tri[tri[, "tr3"] %in% omit.elems, c("tr3", "arc3")] <- 0
      new.tri <- tri[!(tri[,"elem"] %in% omit.elems), ]
    }
    new.tri
  }

  ## Main program

  if (!requireNamespace("tripack", quietly=TRUE))
    stop()

  elem.build <- matrix(c(2, 3, 3, 1, 1, 2), nrow=3, ncol=2, byrow=TRUE,
                       dimnames=list(c("arc1", "arc2", "arc3"),
                                     c("pt1", "pt2")))

  tri <- tripack::triangles(mesh)
  tri <- cbind(elem=seq_len(nrow(tri)), tri)

  itr <- 0
  old.tri <- new.tri <- tri
  while (any(!is.na(new.tri)) & itr < max.itr) {
    itr <- itr + 1
    old.tri <- new.tri
    new.tri <- ModTri(old.tri)
    if (inherits(new.tri, "numeric")) {
      old.tri <- t(as.matrix(new.tri))
      new.tri <- NA
    }
  }

  tri <- old.tri
  n <- nrow(tri)

  tri <- tri[tri[, "tr1"] == 0 | tri[, "tr2"] == 0 | tri[, "tr3"] == 0, ]
  elems <- tri[, "elem"]

  tri[tri[, "tr1"] == 0 & (tri[, "tr2"] != 0 & tri[, "tr3"] != 0), "node1"] <- 0
  tri[tri[, "tr2"] == 0 & (tri[, "tr1"] != 0 & tri[, "tr3"] != 0), "node2"] <- 0
  tri[tri[, "tr3"] == 0 & (tri[, "tr1"] != 0 & tri[, "tr2"] != 0), "node3"] <- 0

  arcs <- c()
  for (i in seq_along(tri[,1])) {
    nodes.in.elem <- array(tri[i, c("node1", "node2", "node3")])
    if (0 %in% nodes.in.elem) {
      arcs <- rbind(arcs, nodes.in.elem[nodes.in.elem != 0])
    } else {
      if (tri[i, "tr1"] == 0)
        arcs <- rbind(arcs, nodes.in.elem[c(2, 3)])
      if (tri[i, "tr2"] == 0)
        arcs <- rbind(arcs, nodes.in.elem[c(1, 3)])
      if (tri[i, "tr3"] == 0)
        arcs <- rbind(arcs, nodes.in.elem[c(1, 2)])
    }
  }

  sort.arcs <- arcs[1, ]
  arcs <- arcs[-1, ]

  for (i in seq_along(arcs[,1])) {
    node <- sort.arcs[length(sort.arcs)]
    has.shared.node <- arcs[, 1] %in% node | arcs[, 2] %in% node

    adjacent.arc <- arcs[has.shared.node, ]

    if (is.na(adjacent.arc[1]))
      return(NULL)

    if (node != adjacent.arc[1])
      adjacent.arc <- rev(adjacent.arc)
    sort.arcs <- c(sort.arcs, adjacent.arc)
    arcs <- arcs[!has.shared.node, ]
    if (inherits(arcs, "numeric"))
      arcs <- t(as.matrix(arcs))
  }

  sort.arcs <- unique(sort.arcs)

  x <- mesh$x[sort.arcs]
  y <- mesh$y[sort.arcs]

  return(as(structure(c(x, y), .Dim=c(length(x), 2)), "gpc.poly"))
}
