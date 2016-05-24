#' @export
FWE <- function(PV, stepdown=TRUE, cfun=NULL) {
  B <- dim(PV)[1] - 1
  P <- dim(PV)[2]
  p.adj <- array(0, dim=c(P, 1))
  if (stepdown) { ## Stepdown shortcut (requires minumum CF)
    cat('\nStepdown Min-P\n')
    p.obs <- PV[1, ] ## raw p-values of the observed statistics
    p.ord <- sort(p.obs, decreasing=FALSE) ## sorted in ascending order
    o <- order(p.obs, decreasing=FALSE)
    PV.ord <- PV[, o]
    ts <- apply(PV.ord, 1, function (p) -min(p))
    p.adj[1] <- max(p.ord[1], mean.default(ts[-1] >= ts[1]))
    if (P > 2) {
      for (j in 2:(P - 1)) {
        ts.j <- apply(PV.ord[, j:P], 1, function (p) -min(p))
        p.adj[j] <- max(mean.default(ts.j[-1] >= ts.j[1]),
                        p.adj[j - 1])
        p.adj[j] <- max(p.ord[j], p.adj[j])
      }
    }
    p.adj[P] <- max(p.ord[P], p.adj[P - 1])
    p.adj[o] <- p.adj
  } else { ## Closed testing (enumerate all intersections - slow!)
    cat('\nClosed Testing\n')
    resamp <- dim(PV)[1] ## number of permutations(?)
    k <- dim(PV)[2] ## number of tests
    rows <- 2^k
    ncycles <- rows
    x <- matrix(0, rows, k) ## indicates tests in each intersection
    for (i in 1:k) { ## loop over tests
      ncycles <- ncycles/2
      nreps <- rows/(2*ncycles)
      zo <- matrix(c(0, 1), nreps, 2, byrow=TRUE)
      zoc <- rbind(as.matrix(zo[, 1]), as.matrix(zo[, 2]))
      settings <- matlab::repmat(zoc, c(1, ncycles))
      x[, k - i + 1] <- settings
    }
    x <- x[-1, ]
    PV2 <- matrix(, resamp, (rows - 1))
    for (j in 1:(rows - 1)) { ## computationally intensive
      if (j %% 50 == 1) cat('Intersection', j, 'of', rows - 1, '\n')
      PV2[, j] <- apply(as.matrix(PV[, x[j, ]==1]), 1, cfun, B=B)
    }
    rawP <- apply(PV2, 2, function(c) mean(c[-1] >= c[1], na.rm=TRUE))
    for (l in 1:k) { ## loop over tests
      ## adjusted p-value is maximum of p-values of intersection
      ## hypotheses to which it belongs
      p.adj[l] <- max(rawP[x[, l]==1])
    }
  }
  rownames(p.adj) <- colnames(PV)
  return(p.adj)
}
