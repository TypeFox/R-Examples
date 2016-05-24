# This file contains custom statistics to be used in conjunction with the 
# gof function. Use these examples to build your own statistics for use with 
# the gof function. Moreover, some helper functions for the statistics can be 
# found here, e.g., functions for computing ROC, PR, and AUC measures based 
# on the ROCR package.


# GOF function for computing dyad-wise shared partner statistics
dsp <- function(mat) {
  d <- summary(mat ~ dsp(0:(nrow(mat) - 2)))
  names(d) <- 0:(length(d) - 1)
  attributes(d)$label <- "Dyad-wise shared partners"
  return(d)
}

# GOF function for computing edge-wise shared partner statistics
esp <- function(mat) {
  d <- summary(mat ~ esp(0:(nrow(mat) - 2)))
  names(d) <- 0:(length(d) - 1)
  attributes(d)$label <- "Edge-wise shared partners"
  return(d)
}

# GOF function for computing non-edge-wise shared partner statistics
nsp <- function(mat) {
  d <- summary(mat ~ nsp(0:(nrow(mat) - 2)))
  names(d) <- 0:(length(d) - 1)
  attributes(d)$label <- "Non-edge-wise shared partners"
  return(d)
}

# GOF function for computing the degree distribution
deg <- function(mat) {
  d <- summary(network(as.matrix(mat), directed = FALSE) ~ degree(0:(nrow(mat) 
      - 1)))
  names(d) <- 0:(length(d) - 1)
  attributes(d)$label <- "Degree"
  return(d)
}

# GOF function for computing the degree distribution for the first mode
b1deg <- function(mat) {
  d <- summary(network(as.matrix(mat), directed = FALSE, bipartite = TRUE) ~ 
      b1degree(0:nrow(mat)))
  names(d) <- 0:(length(d)- 1)
  attributes(d)$label <- "Degree (first mode)"
  return(d)
}

# GOF function for computing the degree distribution for the second mode
b2deg <- function(mat) {
  d <- summary(network(as.matrix(mat), directed = FALSE, bipartite = TRUE) ~ 
      b2degree(0:ncol(mat)))
  names(d) <- 0:(length(d)- 1)
  attributes(d)$label <- "Degree (second mode)"
  return(d)
}

# GOF function for computing the degree distribution
odeg <- function(mat) {
  d <- summary(network(as.matrix(mat), directed = TRUE) ~ odegree(0:(nrow(mat) 
      - 1)))
  names(d) <- 0:(length(d) - 1)
  attributes(d)$label <- "Outdegree"
  return(d)
}

# GOF function for computing the degree distribution
ideg <- function(mat) {
  d <- summary(network(as.matrix(mat), directed = TRUE) ~ idegree(0:(nrow(mat) 
      - 1)))
  names(d) <- 0:(length(d) - 1)
  attributes(d)$label <- "Indegree"
  return(d)
}

# GOF function for computing the k-star distribution
kstar <- function(mat) {
  d <- summary(network(as.matrix(mat), directed = FALSE) ~ kstar(0:(nrow(mat) 
      - 1)))
  names(d) <- 0:(length(d) - 1)
  attributes(d)$label <- "k-star"
  return(d)
}

# GOF function for computing the k-star distribution on the first mode
b1star <- function(mat) {
  d <- summary(network(as.matrix(mat), directed = FALSE, bipartite = TRUE) ~ 
      b1star(0:nrow(mat)))
  names(d) <- 0:(length(d) - 1)
  attributes(d)$label <- "k-star (first mode)"
  return(d)
}

# GOF function for computing the k-star distribution on the second mode
b2star <- function(mat) {
  d <- summary(network(as.matrix(mat), directed = FALSE, bipartite = TRUE) ~ 
      b2star(0:nrow(mat)))
  names(d) <- 0:(length(d) - 1)
  attributes(d)$label <- "k-star (second mode)"
  return(d)
}

# GOF function for computing the outgoing k-star distribution
ostar <- function(mat) {
  d <- summary(network(as.matrix(mat), directed = TRUE) ~ ostar(0:(nrow(mat) 
      - 1)))
  names(d) <- 0:(length(d) - 1)
  attributes(d)$label <- "Outgoing k-star"
  return(d)
}

# GOF function for computing the incoming k-star distribution
istar <- function(mat) {
  d <- summary(network(as.matrix(mat), directed = TRUE) ~ istar(0:(nrow(mat) 
      - 1)))
  names(d) <- 0:(length(d) - 1)
  attributes(d)$label <- "Incoming k-star"
  return(d)
}

# GOF function for computing the degree distribution
kcycle <- function(mat) {
  d <- summary(mat ~ cycle(0:(nrow(mat) - 1)))
  names(d) <- 0:(length(d) - 1)
  attributes(d)$label <- "Cycle"
  return(d)
}

# GOF function for computing geodesic distances
geodesic <- function(mat) {
  mat[is.na(mat)] <- 0
  fillup <- function(x, another.length) {  # fill up x if shorter
    difference <- length(x) - another.length
    inf.value <- x[length(x)]
    if (difference < 0) {  # x is shorter
      x <- x[1:(length(x) - 1)]
      x <- c(x, rep(0, abs(difference)), inf.value)
    } else if (difference > 0) {
      x <- x[1:(length(x) - difference)]
      x <- c(x, inf.value)
    }
    return(x)
  }
  g <- fillup(ergm.geodistdist(network(as.matrix(mat), directed = TRUE)), 
      nrow(mat) - 1)
  attributes(g)$label <- "Geodesic distances"
  return(g)
}

# GOF function for computing triad census statistics in directed graphs
triad.directed <- function(mat) {
  tr <- sna::triad.census(network(as.matrix(mat), directed = TRUE), 
      mode = "digraph")[1, ]
  attributes(tr)$label <- "Triad census"
  return(tr)
}

# GOF function for computing triad census statistics in undirected graphs
triad.undirected <- function(mat) {
  tr <- sna::triad.census(network(as.matrix(mat), directed = FALSE), 
      mode = "graph")[1, ]
  attributes(tr)$label <- "Triad census"
  return(tr)
}

# helper function: create community comembership matrix
comemb <- function(vec) {
  comemb <- matrix(0, nrow = length(vec), ncol = length(vec))
  for (a in 1:length(vec)) {
    for (b in 1:length(vec)) {
      if (vec[a] == vec[b]) {
        comemb[a, b] <- 1
      } else {
        comemb[a, b] <- 0
      }
    }
  }
  return(comemb)
}

# GOF function for computing Walktrap modularity distribution
walktrap.modularity <- function(mat) {
  mat[is.na(mat)] <- 0
  if (xergm.common::is.mat.directed(as.matrix(mat))) {
    m <- "directed"
  } else {
    m <- "undirected"
  }
  if (sum(mat) == 0) {
    mod <- 0
  } else {
    g <- igraph::graph.adjacency(as.matrix(mat), mode = m)
    wt <- igraph::walktrap.community(g)
    mod <- igraph::modularity(wt)
  }
  attributes(mod)$label <- "Modularity (walktrap)"
  return(mod)
}

# ROC for Walktrap community detection algorithm
walktrap.roc <- function(sim, obs) {
  for (i in 1:length(sim)) {
    sim[[i]][is.na(sim[[i]])] <- 0
  }
  for (i in 1:length(obs)) {
    obs[[i]][is.na(obs[[i]])] <- 0
  }
  fun <- function(x) {
    m <- xergm.common::is.mat.directed(as.matrix(x))
    if (m == TRUE) {
      m <- "directed"
    } else {
      m <- "undirected"
    }
    if (sum(x) == 0) {
      memb <- rep(1, nrow(x))
    } else {
      g <- igraph::graph.adjacency(as.matrix(x), mode = m)
      memb <- igraph::walktrap.community(g)$membership
    }
    return(comemb(memb))
  }
  sim <- lapply(sim, fun)
  obs <- lapply(obs, fun)
  object <- suppressMessages(rocpr(sim, obs)[-c(4, 5, 8, 9)])
  class(object) <- "roc"
  object$label <- "Walktrap community comembership prediction"
  attributes(object)$label <- object$label
  return(object)
}

# PR for Walktrap community detection algorithm
walktrap.pr <- function(sim, obs) {
  for (i in 1:length(sim)) {
    sim[[i]][is.na(sim[[i]])] <- 0
  }
  for (i in 1:length(obs)) {
    obs[[i]][is.na(obs[[i]])] <- 0
  }
  fun <- function(x) {
    m <- xergm.common::is.mat.directed(as.matrix(x))
    if (m == TRUE) {
      m <- "directed"
    } else {
      m <- "undirected"
    }
    if (sum(x) == 0) {
      memb <- rep(1, nrow(x))
    } else {
      g.obs <- igraph::graph.adjacency(as.matrix(x), mode = m)
      memb <- igraph::walktrap.community(g.obs)$membership
    }
    return(comemb(memb))
  }
  sim <- lapply(sim, fun)
  obs <- lapply(obs, fun)
  object <- suppressMessages(rocpr(sim, obs)[-c(2, 3, 6, 7)])
  class(object) <- "pr"
  object$label <- "Walktrap community comembership prediction"
  attributes(object)$label <- object$label
  return(object)
}


# GOF function for computing fast and greedy modularity distribution
fastgreedy.modularity <- function(mat) {
  mat[is.na(mat)] <- 0
  if (sum(mat) == 0) {
    mod <- 0
  } else {
    g <- igraph::graph.adjacency(mat, mode = "undirected")
    wt <- igraph::fastgreedy.community(g)
    mod <- igraph::modularity(wt)
  }
  attributes(mod)$label <- "Modularity (fast & greedy)"
  return(mod)
}

# ROC for fast & greedy community detection algorithm
fastgreedy.roc <- function(sim, obs) {
  for (i in 1:length(sim)) {
    sim[[i]][is.na(sim[[i]])] <- 0
  }
  for (i in 1:length(obs)) {
    obs[[i]][is.na(obs[[i]])] <- 0
  }
  fun <- function(x) {
    if (sum(x) == 0) {
      memb <- rep(1, nrow(x))
    } else {
      g <- igraph::graph.adjacency(as.matrix(x), mode = "undirected")
      memb <- igraph::fastgreedy.community(g)$membership
    }
    return(comemb(memb))
  }
  sim <- lapply(sim, fun)
  obs <- lapply(obs, fun)
  object <- suppressMessages(rocpr(sim, obs)[-c(4, 5, 8, 9)])
  class(object) <- "roc"
  object$label <- "Fast $ greedy community comembership prediction"
  attributes(object)$label <- object$label
  return(object)
}

# PR for fast & greedy community detection algorithm
fastgreedy.pr <- function(sim, obs) {
  for (i in 1:length(sim)) {
    sim[[i]][is.na(sim[[i]])] <- 0
  }
  for (i in 1:length(obs)) {
    obs[[i]][is.na(obs[[i]])] <- 0
  }
  fun <- function(x) {
    if (sum(x) == 0) {
      memb <- rep(1, nrow(x))
    } else {
      g <- igraph::graph.adjacency(as.matrix(x), mode = "undirected")
      memb <- igraph::fastgreedy.community(g)$membership
    }
    return(comemb(memb))
  }
  sim <- lapply(sim, fun)
  obs <- lapply(obs, fun)
  object <- suppressMessages(rocpr(sim, obs)[-c(2, 3, 6, 7)])
  class(object) <- "pr"
  object$label <- "Fast $ greedy community comembership prediction"
  attributes(object)$label <- object$label
  return(object)
}

# GOF function for computing maximal modularity distribution
maxmod.modularity <- function(mat) {
  mat[is.na(mat)] <- 0
  if (xergm.common::is.mat.directed(as.matrix(mat))) {
    m <- "directed"
  } else {
    m <- "undirected"
  }
  if (sum(mat) == 0) {
    mod <- 0
  } else {
    g <- igraph::graph.adjacency(mat, mode = m)
    wt <- igraph::optimal.community(g)
    mod <- igraph::modularity(wt)
  }
  attributes(mod)$label <- "Maximum modularity"
  return(mod)
}

# ROC for maximal modularity community detection algorithm
maxmod.roc <- function(sim, obs) {
  for (i in 1:length(sim)) {
    sim[[i]][is.na(sim[[i]])] <- 0
  }
  for (i in 1:length(obs)) {
    obs[[i]][is.na(obs[[i]])] <- 0
  }
  fun <- function(x) {
    m <- xergm.common::is.mat.directed(as.matrix(x))
    if (m == TRUE) {
      m <- "directed"
    } else {
      m <- "undirected"
    }
    if (sum(x) == 0) {
      memb <- rep(1, nrow(x))
    } else {
      g <- igraph::graph.adjacency(as.matrix(x), mode = m)
      memb <- igraph::optimal.community(g)$membership
    }
    return(comemb(memb))
  }
  sim <- lapply(sim, fun)
  obs <- lapply(obs, fun)
  object <- suppressMessages(rocpr(sim, obs)[-c(4, 5, 8, 9)])
  class(object) <- "roc"
  object$label <- "Maximum modularity community comembership prediction"
  attributes(object)$label <- object$label
  return(object)
}

# PR for maximal modularity community detection algorithm
maxmod.pr <- function(sim, obs) {
  for (i in 1:length(sim)) {
    sim[[i]][is.na(sim[[i]])] <- 0
  }
  for (i in 1:length(obs)) {
    obs[[i]][is.na(obs[[i]])] <- 0
  }
  fun <- function(x) {
    m <- xergm.common::is.mat.directed(as.matrix(x))
    if (m == TRUE) {
      m <- "directed"
    } else {
      m <- "undirected"
    }
    if (sum(x) == 0) {
      memb <- rep(1, nrow(x))
    } else {
      g <- igraph::graph.adjacency(as.matrix(x), mode = m)
      memb <- igraph::optimal.community(g)$membership
    }
    return(comemb(memb))
  }
  sim <- lapply(sim, fun)
  obs <- lapply(obs, fun)
  object <- suppressMessages(rocpr(sim, obs)[-c(2, 3, 6, 7)])
  class(object) <- "pr"
  object$label <- "Maximum modularity community comembership prediction"
  attributes(object)$label <- object$label
  return(object)
}

# GOF function for computing edge betweenness modularity distribution
edgebetweenness.modularity <- function(mat) {
  mat[is.na(mat)] <- 0
  if (xergm.common::is.mat.directed(as.matrix(mat))) {
    m <- "directed"
  } else {
    m <- "undirected"
  }
  if (sum(mat) == 0) {
    mod <- 0
  } else {
    g <- igraph::graph.adjacency(mat, mode = m)
    eb <- igraph::edge.betweenness.community(g)
    mod <- igraph::modularity(eb)
  }
  attributes(mod)$label <- "Modularity (edge betweenness)"
  return(mod)
}

# ROC for edge betweenness community detection algorithm
edgebetweenness.roc <- function(sim, obs) {
  for (i in 1:length(sim)) {
    sim[[i]][is.na(sim[[i]])] <- 0
  }
  for (i in 1:length(obs)) {
    obs[[i]][is.na(obs[[i]])] <- 0
  }
  fun <- function(x) {
    m <- xergm.common::is.mat.directed(as.matrix(x))
    if (m == TRUE) {
      m <- "directed"
    } else {
      m <- "undirected"
    }
    if (sum(x) == 0) {
      memb <- rep(1, nrow(x))
    } else {
      g <- igraph::graph.adjacency(as.matrix(x), mode = m)
      memb <- igraph::edge.betweenness.community(g)$membership
    }
    return(comemb(memb))
  }
  sim <- lapply(sim, fun)
  obs <- lapply(obs, fun)
  object <- suppressMessages(rocpr(sim, obs)[-c(4, 5, 8, 9)])
  class(object) <- "roc"
  object$label <- "Edge betweenness community comembership prediction"
  attributes(object)$label <- object$label
  return(object)
}

# PR for edge betweenness community detection algorithm
edgebetweenness.pr <- function(sim, obs) {
  for (i in 1:length(sim)) {
    sim[[i]][is.na(sim[[i]])] <- 0
  }
  for (i in 1:length(obs)) {
    obs[[i]][is.na(obs[[i]])] <- 0
  }
  fun <- function(x) {
    m <- xergm.common::is.mat.directed(as.matrix(x))
    if (m == TRUE) {
      m <- "directed"
    } else {
      m <- "undirected"
    }
    if (sum(x) == 0) {
      memb <- rep(1, nrow(x))
    } else {
      g <- igraph::graph.adjacency(as.matrix(x), mode = m)
      memb <- igraph::edge.betweenness.community(g)$membership
    }
    return(comemb(memb))
  }
  sim <- lapply(sim, fun)
  obs <- lapply(obs, fun)
  object <- suppressMessages(rocpr(sim, obs)[-c(2, 3, 6, 7)])
  class(object) <- "pr"
  object$label <- "Edge betweenness community comembership prediction"
  attributes(object)$label <- object$label
  return(object)
}

# GOF function for computing spinglass modularity distribution
spinglass.modularity <- function(mat) {
  mat[is.na(mat)] <- 0
  if (xergm.common::is.mat.directed(as.matrix(mat))) {
    m <- "directed"
  } else {
    m <- "undirected"
  }
  if (sum(mat) == 0) {
    mod <- 0
  } else {
    g <- igraph::graph.adjacency(mat, mode = m)
    eb <- igraph::spinglass.community(g)
    mod <- igraph::modularity(eb)
  }
  attributes(mod)$label <- "Modularity (spinglass)"
  return(mod)
}

# ROC for spinglass community detection algorithm
spinglass.roc <- function(sim, obs) {
  for (i in 1:length(sim)) {
    sim[[i]][is.na(sim[[i]])] <- 0
  }
  for (i in 1:length(obs)) {
    obs[[i]][is.na(obs[[i]])] <- 0
  }
  fun <- function(x) {
    m <- xergm.common::is.mat.directed(as.matrix(x))
    if (m == TRUE) {
      m <- "directed"
    } else {
      m <- "undirected"
    }
    if (sum(x) == 0) {
      memb <- rep(1, nrow(x))
    } else {
      g <- igraph::graph.adjacency(as.matrix(x), mode = m)
      memb <- igraph::spinglass.community(g)$membership
    }
    return(comemb(memb))
  }
  sim <- lapply(sim, fun)
  obs <- lapply(obs, fun)
  object <- suppressMessages(rocpr(sim, obs)[-c(4, 5, 8, 9)])
  class(object) <- "roc"
  object$label <- "Spinglass community comembership prediction"
  attributes(object)$label <- object$label
  return(object)
}

# PR for spinglass community detection algorithm
spinglass.pr <- function(sim, obs) {
  for (i in 1:length(sim)) {
    sim[[i]][is.na(sim[[i]])] <- 0
  }
  for (i in 1:length(obs)) {
    obs[[i]][is.na(obs[[i]])] <- 0
  }
  fun <- function(x) {
    m <- xergm.common::is.mat.directed(as.matrix(x))
    if (m == TRUE) {
      m <- "directed"
    } else {
      m <- "undirected"
    }
    if (sum(x) == 0) {
      memb <- rep(1, nrow(x))
    } else {
      g <- igraph::graph.adjacency(as.matrix(x), mode = m)
      memb <- igraph::spinglass.community(g)$membership
    }
    return(comemb(memb))
  }
  sim <- lapply(sim, fun)
  obs <- lapply(obs, fun)
  object <- suppressMessages(rocpr(sim, obs)[-c(2, 3, 6, 7)])
  class(object) <- "pr"
  object$label <- "Spinglass community comembership prediction"
  attributes(object)$label <- object$label
  return(object)
}


# AUC-PR function -- modified version of: https://github.com/ipa-tys/ROCR/pull/2
aucpr <- function(pred, precision, recall) {
  falsepos <- pred@fp
  truepos <- pred@tp
  n.positive <- pred@n.pos
  aucvalues <- numeric()
  for (j in 1:length(precision)) {
    fp <- falsepos[[j]]
    tp <- truepos[[j]]
    n.pos <- n.positive[[j]]
    prec <- precision[[j]]
    rec <- recall[[j]]
    
    # if two points are too distant from each other, we need to
    # correctly interpolate between them. This is done according to
    # Davis & Goadrich,
    #"The Relationship Between Precision-Recall and ROC Curves", ICML'06
    for (i in seq_along(rec[-length(rec)])) {
      if (tp[i + 1] - tp[i] > 2) {
        skew = (fp[i + 1] - fp[i]) / (tp[i + 1] - tp[i])
        x = seq(1, tp[i + 1] - tp[i], by = 1)
        rec <- append(rec, (x + tp[i]) / n.pos, after = i)
        prec <- append(prec, (x + tp[i]) / (tp[i] + fp[i] + x + skew * x), 
            after = i)
      }
    }
    
    auc <- 0
    for (i in 2:length(rec)) {
        auc <- auc + 0.5 * (rec[i] - rec[i-1]) * (prec[i] + prec[i-1])
    }
    aucvalues <- c(aucvalues, auc)
  }
  return(aucvalues)
}

# GOF function for ROC curves and PR curves
rocprgof <- function(sim, obs, pr.impute = "poly4") {
  
  directed <- sapply(obs, xergm.common::is.mat.directed)
  twomode <- !sapply(obs, xergm.common::is.mat.onemode)
  
  # create random graphs with corresponding tie probability of each time step
  nsim <- length(sim) / length(obs)
  rg <- list()
  for (i in 1:length(obs)) {
    rn <- nrow(obs[[i]])
    cn <- ncol(obs[[i]])
    n <- rn * cn
    dens <- sna::gden(network(as.matrix(obs[[i]]), directed = directed[i], 
        bipartite = twomode[i]))
    rlist <- list()
    for (j in 1:nsim) {
      r <- matrix(rbinom(n = n, size = 1, prob = dens), nrow = rn, ncol = cn)
      if (twomode[[i]] == FALSE) {
        diag(r) <- 0
        if (directed[[i]] == FALSE) {
          r <- symmetrize(r, rule = "upper")
        }
      }
      rg[length(rg) + 1] <- Matrix(r)
    }
  }
  
  # ROCR
  target.pr <- list()
  rgraph.pr <- list()
  target.y <- list()
  for (j in 1:length(obs)) {
    net <- obs[[j]]
    index.start <- j * nsim - nsim + 1 # sim 1-100 for obs 1, 101-200 for 2 etc.
    index.stop <- (j + 1) * nsim - nsim
    sums <- 0
    rg.sums <- 0
    for (i in index.start:index.stop) {
      sums <- sums + sim[[i]]
      rg.sums <- rg.sums + rg[[i]]
    }
    sums <- sums / length(sim)
    rg.sums <- rg.sums / length(rg)
    if (directed[[j]] == TRUE || twomode[[j]] == TRUE) {
      pr <- c(as.matrix(sums))
      y <- c(as.matrix(net))
      rg.pr <- c(as.matrix(rg.sums))
    } else {
      pr <- sums[lower.tri(sums)]
      y <- as.matrix(net)[lower.tri(as.matrix(net))]
      rg.pr <- rg.sums[lower.tri(rg.sums)]
    }
    pr <- pr[!is.na(y)]
    rg.pr <- rg.pr[!is.na(y)]
    y <- y[!is.na(y)]
    target.pr[[j]] <- pr
    rgraph.pr[[j]] <- rg.pr
    target.y[[j]] <- y
  }
  pred <- prediction(target.pr, target.y)
  roc <- performance(pred, "tpr", "fpr")  # ROC curve
  pr <- performance(pred, "ppv", "tpr")  # precision-recall curve
  rg.pred <- prediction(rgraph.pr, target.y)
  rg.roc <- performance(rg.pred, "tpr", "fpr")  # ROC curve
  rg.pr <- performance(rg.pred, "ppv", "tpr")  # precision-recall curve
  
  # impute the first PR value (which is sometimes NaN)
  for (j in 1:length(pr@y.values)) {
    fp <- pred@fp[[j]]
    tp <- pred@tp[[j]]
    if (fp[1] == 0 & tp[1] == 0) {
      pr@y.values[[j]][1] <- 1
    } else if (pr.impute == "no") {
      message(paste0("t = ", j, ": warning -- the first PR value was not ", 
        "imputed; this may lead to underestimated AUC-PR values."))
      # do nothing
    } else if (is.nan(pr@y.values[[j]][1])) {
      if (pr.impute == "second") {
        message(paste0("t = ", j, ": imputing the first PR value by the ", 
            "next (= adjacent) value."))
        pr@y.values[[j]][1] <- pr@y.values[[j]][2]
      } else if (pr.impute == "one") {
        message(paste0("t = ", j, ": imputing the first PR value by the ", 
            "maximum value of 1."))
        pr@y.values[[j]][1] <- 1
      } else if (grepl("^poly[1-9]", pr.impute)) {
        num <- as.numeric(substr(pr.impute, 5, 5))
        message(paste0("t = ", j, ": imputing the first PR value using a ", 
            "polynomial of order ", num, ". Check the results by plotting ",
            "the GOF object using the \"pr.poly = ", num, "\" argument."))
        p <- data.frame(poly(pr@x.values[[j]], num, raw = TRUE))
        fit <- lm(pr@y.values[[j]] ~ ., data = p)
        pr@y.values[[j]][1] <- predict(fit, newdata = p[1, ])
      } else {
        message(paste0("t = ", j, ": PR imputation method not recognized. ", 
            "Not using any imputation."))
      }
      if (pr@y.values[[j]][1] < 0) {
        pr@y.values[[j]][0] <- 0
      }
      if (pr@y.values[[j]][1] > 1) {
        pr@y.values[[j]][1] <- 1
      }
    }
  }
  
  auc.roc <- unlist(performance(pred, measure = "auc")@y.values)  # ROC-AUC
  auc.pr <- aucpr(pred, precision = pr@y.values, recall = pr@x.values)  # PR-AUC
  rg.auc.roc <- unlist(performance(rg.pred, measure = "auc")@y.values)
  rg.auc.pr <- aucpr(rg.pred, precision = rg.pr@y.values, 
      recall = rg.pr@x.values)
  
  rocpr <- list()
  rocpr$type <- "rocpr"
  rocpr$auc.roc <- auc.roc
  rocpr$auc.roc.rgraph <- rg.auc.roc
  rocpr$auc.pr <- auc.pr
  rocpr$auc.pr.rgraph <- rg.auc.pr
  rocpr$roc <- roc
  rocpr$roc.rgraph <- rg.roc
  rocpr$pr <- pr
  rocpr$pr.rgraph <- rg.pr
  class(rocpr) <- "rocpr"
  return(rocpr)
}

# wrapper function for rocprgof without pr.impute argument (for use with gof)
rocpr <- function(sim, obs) {
  object <- suppressMessages(rocprgof(sim, obs))
  object$label <- "ROC and PR curve"
  attributes(object)$label <- object$label
  return(object)
}

# wrapper function for ROC only (for use with gof)
roc <- function(sim, obs) {
  object <- suppressMessages(rocpr(sim, obs)[-c(4, 5, 8, 9)])
  class(object) <- "roc"
  object$label <- "Receiver-operating characteristics"
  attributes(object)$label <- object$label
  return(object)
}

# wrapper function for PR only (for use with gof)
pr <- function(sim, obs) {
  object <- suppressMessages(rocpr(sim, obs)[-c(2, 3, 6, 7)])
  class(object) <- "pr"
  object$label <- "Precision-recall curve"
  attributes(object)$label <- object$label
  return(object)
}
