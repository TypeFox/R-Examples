# Total Interaction Indices (TII) by the pick-freeze method 
# Author: Jana Fruth (2014)

sobolTIIpf <- function (model = NULL, X1, X2, ...){
  if ((ncol(X1) != ncol(X2)) | (nrow(X1) != nrow(X2))) {
    stop("The samples X1 and X2 must have the same dimensions")
  }
  p <- ncol(X1)
  indices.list <- subsets(set = 1:p, size = 1)
  X <- X1
  # X is the final matrix to be evaluated
  # contains: X1, all X1 with col i from X2
  for (i in indices.list) {
    Xb <- X1
    Xb[, i] <- X2[, i]
    X <- rbind(X, Xb)
  }
  x <- list(model = model, X1 = X1, X2 = X2, X = X, call = match.call())
  # match.call erinnert sich nur an den ursprunglichen Funktionsaufruf (superliuowen.. usw)
  class(x) <- "sobolTIIpf"
  if (!is.null(x$model)) {
    response(x, ...)
    # response rechnet das Modell an den ganzen x aus und hangt es als y dran
    tell(x, ...)
  }
  return(x)
}



tell.sobolTIIpf <- function (x, y = NULL, ...) {
  id <- deparse(substitute(x))
  if (!is.null(y)) {
    x$y <- y
  }
  else if (is.null(x$y)) {
    stop("y not found")
  }
  p <- ncol(x$X1)
  n <- nrow(x$X1)
  indices.list2 <- subsets(set = 1:p, size = 2)[-(1:p)]
  ni <- length(indices.list2)
  indices.labels <- lapply(indices.list2, function(i) paste(colnames(x$X1)[i], 
                                                            collapse = "*"))
  data <- matrix(x$y, nrow = n)
  # y vector seperated columnwise according to the matrices in X
  V <- data.frame(original = estim.sobolTIIpf(data, indices.list2 = indices.list2))
  tii <- V[2:(ni+1),,drop=FALSE]
  colnames(tii) <- c("original")
  rownames(tii) <- indices.labels
  x$V <- V[1,1]
  x$tii.unscaled <- tii
  x$tii.scaled <- tii[,1,drop=FALSE]/x$V
  assign(id, x, parent.frame())
}


estim.sobolTIIpf <- function (data, i = 1:nrow(data), indices.list2) 
{
  d <- as.matrix(data[i, ])
  n <- nrow(d)
  ni <- length(indices.list2)
  p <- indices.list2[[ni]][2]
  V <- var(as.numeric(data)) * (n - 1)/n
  tii <- numeric(ni)  
  DC.but.i <- numeric(p)
  for (i in 1:p) {
    mu <- 1/2 * mean(data[,1] + data[, i+1])
    DC.but.i[i] <- mean(data[, i+1] * data[,1]) - mu^2
  }
  DC.but.ij <- numeric(ni)
  for (r in 1:ni) {
    i <- indices.list2[[r]][1]
    j <- indices.list2[[r]][2]
    mu <- 1/2 * mean(data[, i+1] + data[, j+1])
    DC.but.ij[r] <- mean(data[, i+1] * data[, j+1]) - mu^2
  }
  for (r in 1:ni) {
    tii[r] <- V + DC.but.ij[r] - sum(DC.but.i[indices.list2[[r]]])
  }
  return(c(V, tii))
}

print.sobolTIIpf <- function (x, ...) 
{
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if (!is.null(x$y)) {
    cat("\nModel runs:", length(x$y), "\n")
    if (!is.null(x$tii.scaled)) {
      cat("\nscaled tii\n")
      print(x$tii.scaled)
    }
  }
  else {
    cat("(empty)\n")
  }
}

plot.sobolTIIpf <- function(x, ylim = NULL, ...)
{
  if (is.null(ylim)){
    ylim <- range(x$tii.scaled[,1])
  }
  if (!is.null(x$y)) {
    nodeplot(x$tii.scaled, ylim = ylim)
  }
}

plotFG.sobolTIIpf <- function (x) 
{
  if (!is.null(x$y)) {
    max.thickness <- 15
    diameter <- 28
    p <- ncol(x$X1)
    tii <- x$tii.unscaled[,1]
    active <- which(tii > 0)
    tii <- tii[active]
    E <- t(combn(p, 2))
    E <- E[active, ]
    if (requireNamespace("igraph", quietly = TRUE)){ 
      g <- igraph::graph(as.vector(t(E)), n = p, directed = FALSE)
    }
    if (requireNamespace("igraph", quietly = TRUE)){ 
      layout <- igraph::layout.fruchterman.reingold(g)
    }
    names <- colnames(x$X1)
    edge.weight.scale <- tii/(max(tii)) * max.thickness
    if (requireNamespace("igraph", quietly = TRUE)){ 
      igraph::plot.igraph(g,layout = layout, edge.width = edge.weight.scale, vertex.frame.color="darkgrey",
                          vertex.color = "white", vertex.label = names, vertex.size=diameter)
    }
  }
}
