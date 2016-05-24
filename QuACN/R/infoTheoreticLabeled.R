infoTheoreticLabeledV1 <- function(g, dist=NULL, coeff="lin", custCoeff=NULL, coeffMatrix=NULL, lambda=1000) {
  if (class(g)[1] != "graphNEL")
    stop("'g' has to be a 'graphNEL' object")
  stopifnot(.validateGraph(g))
  
  if (is.null(dist))
    dist <- distanceMatrix(g)

  labels <- .nodeDataVector(g, "atom->symbol")
  uniq.labels <- unique(labels)

  # if there is no coefficient matrix, create one from coeff/custCoeff
  # and the masses of the atoms in the graph
  if (is.null(coeffMatrix)) {
    ##data(sysdata, envir=environment())
    diam <- max(dist)
    ci <- .infoTheoreticCoeff(coeff, custCoeff, diam)
    coeffMatrix <- matrix(0, nrow=diam, ncol=length(uniq.labels))
    colnames(coeffMatrix) <- uniq.labels
    uranium <- round(.chemicalElements()[92, "mass"])
    relmasses <- sapply(uniq.labels, function(label) {
      absmass <- .chemicalElements()[.chemicalElements()$symbol == label, "mass"]
      round(absmass) / uranium
    })

    for (i in 1:diam)
      for (label in uniq.labels)
        coeffMatrix[i, label] <- ci[[i]] - relmasses[[label]]
  }

  # calculate the functional by looking at each pair of vertices
  # and weighting them according to distance (which j-sphere it is on) and label
  fvi <- sapply(nodes(g), function(vi) {
    sum(sapply(nodes(g), function(target) {
      if (target == vi)
        0
      else {
        k <- dist[vi, target]
        l <- labels[[target]]
        coeffMatrix[k, l]
      }
    }))
  })

  .infoTheoretic(fvi, lambda)
}

infoTheoreticLabeledV2 <- function(g, ci=NULL, lambda=1000) {
  if (class(g)[1] != "graphNEL")
    stop("'g' has to be a 'graphNEL' object")

  labels <- .nodeDataVector(g, "atom->symbol")
  uniq.labels <- unique(labels)
  if (is.null(ci)) {
    ci <- rep(1, length(uniq.labels))
    names(ci) <- uniq.labels
  }
  else if (is.list(ci)) {
    tmp <- as.numeric(ci)
    names(tmp) <- names(ci)
    ci <- tmp
  }
  ci <- ci[uniq.labels]

  ig <- .G2IG(g)
  n <- numNodes(g)

  # determine number of all possible shortest paths
  fvi <- sapply(1:n, function(vi) {
    asp <- igraph::get.all.shortest.paths(ig, from=vi)$res
    asp.lns <- sapply(asp, length)

    # which nodes are in local information graphs of length ll?
    nod.loc <- sapply(2:max(asp.lns), function(ll) {
      unique(unlist(asp[asp.lns == ll]))
    }, simplify=FALSE)

    # match them with the labels
    nod.loc.lab <- sapply(nod.loc, function(nl) {
      a <- sort(as.numeric(nl))
      labels[a]
    }, simplify=FALSE)

    # count labels for local information graphs of each length
    fvi.detail <- sapply(nod.loc.lab, function(nll) {
      sapply(uniq.labels,function(tl) {
	length(nll[nll == tl])
      })
    }, simplify=FALSE)
    fvi.detail <- do.call(cbind, fvi.detail)

    # sum it up for each label
    sums <- apply(fvi.detail, 1, sum)

    sum(ci*sums[match(names(ci), names(sums))])
  })

  names(fvi) <- nodes(g)

  .infoTheoretic(fvi, lambda)
}

infoTheoreticLabeledE <- function(g, dist=NULL, coeff="lin", custCoeff=NULL, lambda=1000) {
  if (class(g)[1] != "graphNEL")
    stop("'g' has to be a 'graphNEL' object")
  if (is.null(dist))
    dist <- distanceMatrix(g)

  ci <- .infoTheoreticCoeff(coeff, custCoeff, max(dist))

  bonds <- .edgeDataMatrix(g, "bond")

  ig <- .G2IG(g)
  nam <- nodes(g)
  n <- length(nam)

  fvi <- sapply(1:n, function(vi) {
    # create local information graphs
    asp <- igraph::get.all.shortest.paths(ig, from=vi)$res
    asp.lns <- sapply(asp, length)
    asp <- asp[order(asp.lns)]
    asp.lns <- sapply(asp, length)
    asp <- asp[asp.lns != 1]

    # get graphNEL names
    loc.nods <- sapply(asp, function(a) nam[a], simplify = FALSE)

    # get edges and replace them with their edge labels
    loc.weights <- sapply(loc.nods, function(ln) {
      sapply(1:(length(ln) - 1), function(i) {
        from <- ln[[i]]
        to <- ln[[i + 1]]
        bonds[from, to]
      })
    })
    loc.weights.len <- sapply(loc.weights, length)

    # sum paths weighted by length
    sum(sapply(sort(unique(loc.weights.len)), function(t) {
      ci[[t]] * sum(unlist(loc.weights[loc.weights.len==t]))
    }))
  })

  names(fvi) <- nam

  .infoTheoretic(fvi, lambda)
}
