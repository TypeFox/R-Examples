library(bpkde)

if(require(KernSmooth)) {
  datasets <- data(package = "bpkde")$results[, 3]

  for(d in datasets) {
    data(list = d)
    mat <- get(d)
    x <- seq(min(mat[, 1]) - 1, max(mat[, 1]) + 1, length = 100)
    y <- seq(min(mat[, 2]) - 1, max(mat[, 2]) + 1, length = 100)
    ks <- KernSmooth:::linbin2D(mat, x, y)
    me <- bpkde:::linbin2D(mat, x, y)
    print(d)
    stopifnot(all.equal(ks, me))
  }
}


