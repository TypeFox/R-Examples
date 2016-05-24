# RG W more memeory efficient
`rg_w` <-
function (nodes = 100, arcs = 300, weights = 1, directed = TRUE, seed = NULL) {
  if(!is.null(seed)) 
    set.seed(as.integer(seed))
  if (arcs > 1) {
    rg_w <- data.frame(
      i = sample(1:nodes, (arcs * 1.5), replace = TRUE), 
      j = sample(1:nodes, (arcs * 1.5), replace = TRUE), 
      w = sample(weights, (arcs * 1.5), replace = TRUE))
    rg_w <- rg_w[rg_w[, 1] != rg_w[, 2], ]
    if (directed) {
      rg_w <- rg_w[!duplicated(rg_w[, 1:2]), ]
      rg_w <- rg_w[sample(1:nrow(rg_w), arcs), ]
    } else {
      rg_w <- rg_w[rg_w[, 1] < rg_w[, 2], ]
      rg_w <- rg_w[!duplicated(rg_w[, 1:2]), ]
      rg_w <- rg_w[sample(1:nrow(rg_w), arcs * 0.5), ]
      rg_w <- rbind(as.matrix(rg_w), cbind(rg_w[, 2], rg_w[, 1], rg_w[, 3]))
    }
  } else {
    rg_w <- stats::runif(nodes^2)
    rg_w <- rg_w < arcs
    rg_w <- matrix(data = rg_w, nrow = nodes, ncol = nodes)
    rg_w <- which(rg_w, arr.ind = TRUE)
    rg_w <- rg_w[rg_w[, 1] != rg_w[, 2], ]
    rg_w <- cbind(rg_w, sample(weights, nrow(rg_w), replace = TRUE))
    if (!directed) {
      rg_w <- rg_w[rg_w[, 1] < rg_w[, 2], ]
      rg_w <- rbind(rg_w, cbind(rg_w[, 2], rg_w[, 1], rg_w[, 3]))
    }
  }
  return(as.tnet(rg_w, type = "weighted one-mode tnet"))
}