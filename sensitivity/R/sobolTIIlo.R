# Total Interaction Indices (TII) by the method of Liu & Owen 
# Author: Jana Fruth (2014)

sobolTIIlo <- function (model = NULL, X1, X2, conf = 0.95, ...){
  if ((ncol(X1) != ncol(X2)) | (nrow(X1) != nrow(X2))) {
    stop("The samples X1 and X2 must have the same dimensions")
  }
  p <- ncol(X1)
  indices.list <- subsets(set = 1:p, size = 2)
  X <- X1
  # X is the final matrix to be evaluated
  # contains: X1, all X1 with col i from X2, all X1 with cols i,j from X2
  for (i in indices.list) {
    Xb <- X1
    Xb[, i] <- X2[, i]
    X <- rbind(X, Xb)
  }
  x <- list(model = model, X1 = X1, X2 = X2, X = X, conf = conf, call = match.call())
  # match.call erinnert sich nur an den ursprunglichen Funktionsaufruf (superliuowen.. usw)
  class(x) <- "sobolTIIlo"
  if (!is.null(x$model)) {
    response(x, ...)
    # response rechnet das Modell an den ganzen x aus und hangt es als y dran
    tell(x, ...)
  }
  return(x)
}



tell.sobolTIIlo <- function (x, y = NULL, ...) {
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
  V <- data.frame(original = estim.sobolTIIlo(data, indices.list2 = indices.list2))
  tii.unscaled <- V[2:(ni+1),,drop=FALSE]
  tii.scaled <- tii.unscaled/V[1,1]
  std.error <- sqrt(V[(ni+2):(2*ni+1),])/sqrt(n)
  lower <- tii.scaled - qnorm(1-(1-x$conf)/2) * std.error
  upper <- tii.scaled + qnorm(1-(1-x$conf)/2) * std.error
  tii.scaled <- cbind(tii.scaled, std.error, lower, upper)
  colnames(tii.scaled) <- c("original", "std.error", "min. c.i.", "max. c.i.")
  rownames(tii.scaled) <- indices.labels
  colnames(tii.unscaled) <- c("original")
  rownames(tii.unscaled) <- indices.labels
  x$V <- V[1,1]
  x$tii.unscaled <- tii.unscaled
  x$tii.scaled <- tii.scaled
  assign(id, x, parent.frame())
}


estim.sobolTIIlo <- function (data, i = 1:nrow(data), indices.list2) 
{
  d <- as.matrix(data[i, ])
  n <- nrow(d)
  ni <- length(indices.list2)
  p <- indices.list2[[ni]][2]
  V <- var(d[, 1])
  totalInt <- numeric(ni)
  varInt <- numeric(ni)
  for (index in 1:ni) {
    i <- indices.list2[[index]][1]
    j <- indices.list2[[index]][2]
    deltasq <- ((data[,p+1+index] - data[,i+1] - data[,j+1] + data[,1]))^2
    totalInt[index] <- mean(deltasq)/4
#    varInt[index] <- var(deltasq)/16 - old variance of unscaled index
    varInt[index] <- var(deltasq/(4*V)) # for the scaled indicex
  }
  return(c(V, totalInt, varInt))
}

print.sobolTIIlo <- function (x, ...) 
{
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if (!is.null(x$y)) {
    cat("\nModel runs:", length(x$y), "\n")
    if (!is.null(x$tii.unscaled)) {
      cat("\nunscaled tii\n")
      print(x$tii.unscaled)
    }
    if (!is.null(x$tii.scaled)) {
      cat("\nscaled tii\n")
      print(x$tii.scaled)
    }
  }
  else {
    cat("(empty)\n")
  }
}


plot.sobolTIIlo <- function(x, ylim = NULL, ...)
{
  if (is.null(ylim)){
    ylim <- c(min(rbind(0,x$tii.scaled["min. c.i."])), max(x$tii.scaled["max. c.i."]))
  }
  if (!is.null(x$y)) {
    nodeplot(x$tii.scaled, ylim = ylim)
  }
}


plotFG.sobolTIIlo <- function (x) 
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