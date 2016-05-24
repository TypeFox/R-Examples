### Version 2.0: codeYs are no longer returned, but codes are a list
### with two elements: X and Y.

"bdk" <- function(data, Y, grid = somgrid(), rlen = 100,
                  alpha = c(0.05, 0.01),
                  radius = quantile(nhbrdist, 0.67) * c(1, -1),
                  xweight = 0.75, contin,
                  toroidal = FALSE, n.hood, keep.data = TRUE)
{
  if (!is.numeric(data))
    stop("Argument data should be numeric")
  data <- as.matrix(data)
  
  nd <- nrow(data)
  nx <- ncol(data)
  
  if (is.factor(Y)) {
    YY <- classvec2classmat(Y)
  } else {
    if (is.vector(Y)) {
      YY <- matrix(Y, ncol=1)
    } else {
      YY <- Y
    }
  }
  ny <- ncol(YY)

  if (nrow(YY) != nd)
    stop("Both data matrices should have the same number of rows")

  if (missing(contin))
    contin <- any(abs(rowSums(YY) - 1) > 1e-8)

  ng <- nrow(grid$pts)
  xdists <- ydists <- rep(0, ng)  
  
  starters <- sample(1:nd, ng, replace = FALSE)
  init <- data[starters, , drop = FALSE]
  codes <- init
  if (!contin) {
    ## rescale to .25 - .75 in order to make class transitions easier
    codeYs <- 0.5 + 0.5*(YY[starters,] - 0.5)
  } else {
    codeYs <- YY[starters,]
  }
  
  if (missing(n.hood)) {
    n.hood <- switch(grid$topo,
                     hexagonal = "circular",
                     rectangular = "square")
  } else {
    n.hood <- match.arg(n.hood, c("circular", "square"))
  }
  grid$n.hood <- n.hood
  nhbrdist <- unit.distances(grid, toroidal)

  if (length(radius) == 1)
    radius <- sort(radius * c(1, -1), decreasing = TRUE)

  changes <- rep(0, rlen*2)

  if (contin) {
    res <- .C("BDK_Eucl",
              data = as.double(data),
              Ys = as.double(YY),
              codes = as.double(codes),
              codeYs = as.double(codeYs),
              nhbrdist = as.double(nhbrdist),
              alpha = as.double(alpha),
              radii = as.double(radius),
              xweight = as.double(xweight),
              changes = as.double(changes),
              xdists = as.double(xdists),
              ydists = as.double(ydists),
              n = as.integer(nd),
              px = as.integer(nx),
              py = as.integer(ny),
              ncodes = as.integer(ng),
              rlen = as.integer(rlen),
              PACKAGE = "kohonen")
  } else {
    res <- .C("BDK_Tani",
              data = as.double(data),
              Ys = as.double(YY),
              codes = as.double(codes),
              codeYs = as.double(codeYs),
              nhbrdist = as.double(nhbrdist),
              alpha = as.double(alpha),
              radius = as.double(radius),
              xweight = as.double(xweight),
              changes = as.double(changes),
              xdists = as.double(xdists),
              ydists = as.double(ydists),
              n = as.integer(nd),
              px = as.integer(nx),
              py = as.integer(ny),
              ncodes = as.integer(ng),
              rlen = as.integer(rlen),
              PACKAGE = "kohonen")
  }
  
  changes <- matrix(res$changes, ncol=2)
  codes <- list(X = matrix(res$codes, nrow(init), ncol(init)),
                Y = matrix(res$codeYs, ng, ny))
  colnames(codes$Y) <- colnames(YY)
  colnames(codes$X) <- colnames(data)

  if (keep.data) {
    mapping <- map.kohonen(list(codes = codes), newdata = data, whatmap = 1)
    
    structure(list(data = data, Y = Y, contin = contin,
                   grid = grid, codes = codes,
                   changes = changes, alpha = alpha,
                   radius = radius, toroidal = toroidal,
                   unit.classif = mapping$unit.classif,
                   distances = mapping$distances, method="bdk"),
              class = "kohonen")
  } else {
    structure(list(contin = contin, grid = grid, codes = codes,
                   changes = changes, alpha = alpha,
                   radius = radius, toroidal = toroidal, method="bdk"),
              class = "kohonen")
  }
}

