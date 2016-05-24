### Compute the height of the rooted tree.
get.rooted.tree.height <- function(rooted.tree, tol = .Machine$double.eps^0.5){
  if(!is.rooted(rooted.tree)){
    stop("A rooted tree is required")
  }
  if(rooted.tree$Nnode > 1){
    if(!is.ultrametric(rooted.tree, tol = tol)){
      stop("A rooted tree is not ultrametric")
    }
  }

  edge.from <- rooted.tree$edge[, 1]
  edge.to <- rooted.tree$edge[, 2]

  edge.id <- which(edge.from == (length(rooted.tree$tip.label) + 1))
  n.from <- edge.to[edge.id[1]]
  ret <- rooted.tree$edge.length[edge.id[1]]
  repeat{
    edge.id <- which(edge.from == n.from)
    if(length(edge.id) == 2){
      n.from <- edge.to[edge.id[1]]
      ret <- ret + rooted.tree$edge.length[edge.id[1]]
    } else{
      break
    }
  }
  ret
} # End of get.rooted.tree.height().


### Rescale the rooted tree branche lengths by a constant.
rescale.rooted.tree <- function(rooted.tree, scale.height = 1){
  if(!is.rooted(rooted.tree)){
    stop("A rooted tree is required")
  }
  if(scale.height <= 0){
    stop("The scale.height is not correct.")
  }

  rooted.tree$edge.length <- rooted.tree$edge.length * scale.height
  rooted.tree
} # End of rescale.rooted.tree().


### Reshap a rooted tree to a star shape.
as.star.tree <- function(rooted.tree, keep.bifurcation = TRUE){
  if(!is.rooted(rooted.tree)){
    stop("A rooted tree is required")
  }

  n.tip <- length(rooted.tree$tip.label)
  tree.height <- get.rooted.tree.height(rooted.tree)
  if(keep.bifurcation){
    rooted.tree$edge.length[rooted.tree$edge[, 2] > n.tip] <- 0
    rooted.tree$edge.length[rooted.tree$edge[, 2] <= n.tip] <- tree.height
  } else{
    rooted.tree$edge <- rooted.tree$edge[rooted.tree$edge[, 2] <= n.tip,]
    rooted.tree$edge[, 1] <- max(rooted.tree$edge[, 2]) + 1
    rooted.tree$Nnode <- 1
    rooted.tree$edge.length <- rep(tree.height, n.tip)
  }
  rooted.tree
} # End of as.star.tree().


### Generate a star tree.
gen.star.tree <- function(N, total.height = 1){
  if(N <= 1){
    stop("N > 1.")
  }

  ms.star <- ms(N, 1, opts = paste("-T", sep = " "))
  tree.star <- read.tree(text = ms.star[3])
  tree.star$tip.label <- paste("s", 1:N, sep = "")
  tree.star <- as.star.tree(tree.star)
  th <- get.rooted.tree.height(tree.star)
  ret <- rescale.rooted.tree(tree.star, total.height / th)

  ret
} # End of gen.star.tree().

### Generate a tree with unit heights and K clusters.
gen.unit.K <- function(K, N.K, rate.anc = 10, rate.dec = 10){
  if(K <= 1){
    stop("K > 1.")
  }
  if(any(N.K <= 1)){
    stop("All N.K > 1.")
  }
  if(length(N.K) != K){
    stop("The length of N.K is not equal to K.")
  }
  if(rate.anc <= 0 || rate.dec <= 0){
    stop("The rate.anc or rate.dec <= 0.")
  }

  ms.anc <- ms(K, 1, opts = paste("-T -G", rate.anc, sep = " "))
  tree.anc <- read.tree(text = ms.anc[3])
  tree.anc$tip.label <- paste("a", 1:K, sep = "")

  tree.dec <- NULL
  tree.org <- tree.anc
  tree.star <- tree.anc
  for(k in 1:K){
    ms.dec <- ms(N.K[k], 1, opts = paste("-T -G", rate.dec, sep = " "))
    tree.dec[[k]] <- read.tree(text = ms.dec[3])
    tree.dec[[k]]$tip.label <- paste("d", k, ".", 1:N.K[k], sep = "")
    tree.org <- bind.tree(tree.org, tree.dec[[k]],
                          where = which(tree.org$tip.label ==
                                        paste("a", k, sep = "")))
    tree.star <- bind.tree(tree.star, as.star.tree(tree.dec[[k]]),
                           where = which(tree.star$tip.label ==
                                         paste("a", k, sep = "")))
  }

  tree.anc.height <- get.rooted.tree.height(tree.anc)
  tree.dec.height <- unlist(lapply(tree.dec, get.rooted.tree.height))
  tree.dec.mean.height <- mean(tree.dec.height)
  total.height <- tree.anc.height + tree.dec.mean.height
  max.height <- tree.anc.height + max(tree.dec.height)

  tree.equal <- rescale.rooted.tree(tree.anc, 1 / total.height)
  tmp.scale <- tree.dec.mean.height / total.height
  for(k in 1:K){
    tree.dec.equal <- rescale.rooted.tree(tree.dec[[k]],
                                          tmp.scale / tree.dec.height[k])
    tree.equal <- bind.tree(tree.equal, tree.dec.equal,
                            where = which(tree.equal$tip.label ==
                                          paste("a", k, sep = "")))
  }

  tree.max <- rescale.rooted.tree(tree.org, 1 / max.height)
  tree.star <- rescale.rooted.tree(tree.star, 1 / max.height)

  ret <- list(K = K, N.K = N.K, rate.anc = rate.anc, rate.dec = rate.dec,
              height.anc = tree.anc.height,
              height.dec = tree.dec.height,
              anc = tree.anc, dec = tree.dec,
              org = tree.org,
              equal = tree.equal,
              max = tree.max,
              star = tree.star)
  ret
} # End of gen.unit.K().


### Generate star ancestor tree binded with equal star descendant trees.
gen.equal.star.anc.dec <- function(K, N.K, rate.f = 0.5){
  if(length(N.K) != K){
    stop("The length of N.K is not equal to K.")
  }
  if(rate.f <= 0 || rate.f >= 1){
    stop("The rate.f <= 0 or >= 1.")
  }
  if(K <= 0 || length(N.K) != K || any(N.K <= 0)){
    stop("The K or N.K is not consistant.")
  }

  ms.anc <- ms(K, 1, opts = paste("-T -G", 1, sep = " "))
  tree.anc <- read.tree(text = ms.anc[3])
  tree.anc$tip.label <- paste("a", 1:K, sep = "")

  tree.dec <- NULL
  for(k in 1:K){
    ms.dec <- ms(N.K[k], 1, opts = paste("-T -G", 1, sep = " "))
    tree.dec[[k]] <- read.tree(text = ms.dec[3])
    tree.dec[[k]]$tip.label <- paste("d", k, ".", 1:N.K[k], sep = "")
  }

  tmp.tree <- as.star.tree(tree.anc)
  tmp.height <- get.rooted.tree.height(tmp.tree)
  tree.equal.star <- rescale.rooted.tree(tmp.tree, rate.f / tmp.height)

  for(k in 1:K){
    tmp.tree <- as.star.tree(tree.dec[[k]])
    tmp.height <- get.rooted.tree.height(tmp.tree)
    tmp.tree <- rescale.rooted.tree(tmp.tree, (1 - rate.f) / tmp.height)
    tree.equal.star <- bind.tree(tree.equal.star, tmp.tree,
                                 where = which(tree.equal.star$tip.label ==
                                               paste("a", k, sep = "")))
  }

  ret <- list(K = K, N.K = N.K, rate.f = rate.f,
              anc = tree.anc, dec = tree.dec,
              equal.star = tree.equal.star)
  ret
} # End of gen.equal.anc.dec.star().

