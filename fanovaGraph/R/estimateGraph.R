estimateGraph <- function(f.mat, d, q = NULL, q.arg = NULL, 
    n.tot = NULL, method = "LiuOwen", n.lo = NULL, n.mc = NULL,
    n.fast = 500, L = NULL, M = 6, n.pf = NULL, n.main = 1000, confint = TRUE, print.loop.index= FALSE, ...) {
  p <- choose(d, 2)
  if (is.null(q)) {
      q <- rep("qunif", d)
  } else if (length(q) == 1) {
      q <- rep(q, d)
  }
  if (is.null(q.arg)) {
      q.arg <- rep(list(list()), d)
  } else if (all(sapply(q.arg, function(x) is.numeric(x) & length(x)==1))) {
      q.arg <- rep(list(q.arg), d)
  }
  # begin: input test
  if (!(all(sapply(q.arg, is.list)) & length(q.arg)==d))
    stop("q.arg must be a list of lists of quantile functions parameters")
  
  Xtest <- matrix(0.5,2,d)
   for (j in (1:d))
     Xtest[, j] <- do.call(rep(q,2)[j],c(list(p = Xtest[, j]),rep(q.arg,2)[[j]]))
   if (length(f.mat(Xtest,...))!=2)
     stop("f.mat must be a vectorized function (matrix input must be possible)")
  
  if (!(method %in% c("LiuOwen", "FixFast", "RBD", "PickFreeze"))) 
      stop("method must be set to 'LiuOwen', 'FixFast' ,'RBD' or 'PickFreeze'")
  # end: input test
  
  if (method == "RBD") {
      if (!is.null(L) & !is.null(n.tot)) {
          warning("L will be omitted since n.tot is specified")
          L <- round(n.tot/(2 * (p + d)) - 6 * d, 0)
      }
      if (is.null(L) & !is.null(n.tot)) 
          L <- round(n.tot/(2 * (p + d)) - 6 * d, 0)
      if (is.null(L) & is.null(n.tot)) 
          stop("either n.tot or L must be specified")
      tii <- estimateGraphRBD(f.mat, d, q, q.arg, L, M, print.loop.index, ...)
  }
  
  
   if (method == "PickFreeze") {
      if (!is.null(n.pf) & !is.null(n.tot)) {
          warning("n.pf will be omitted since n.tot is specified")
          n.pf <- round(n.tot/(d + 1), 0)
      }
      if (is.null(n.pf) & !is.null(n.tot)) 
          n.pf <- round(n.tot/(d + 1), 0)
      if (is.null(n.pf) & is.null(n.tot)) 
          stop("either n.tot or n.pf must be specified")
      tii <- estimateGraphPickFreeze(f.mat, d, q, q.arg, n.pf, print.loop.index, ...)
  }
  
  if (method == "LiuOwen") {
    if (!is.null(n.lo) & !is.null(n.tot)) {
      warning("n.lo will be omitted since n.tot is specified")
      n.lo <- round(n.tot/(p + d + 1), 0)
    }
    if (is.null(n.lo) & !is.null(n.tot)) 
      n.lo <- round(n.tot/(p + d + 1), 0)
    if (is.null(n.lo) & is.null(n.tot)) 
      stop("either n.tot or n.lo must be specified")
    tii <- estimateGraphLiuOwen(f.mat, d, q, q.arg, n.lo, confint, print.loop.index, ...)
  }
      
  if (method == "FixFast") {
      if (!is.null(n.mc) & !is.null(n.tot)) {
          warning("n.mc will be omitted since n.tot is specified")
          n.mc <- round(n.tot/(p * n.fast), 0)
      }
      if (is.null(n.mc) & !is.null(n.tot)) 
          n.mc <- round(n.tot/(p * n.fast), 0)
      if (is.null(n.mc) & is.null(n.tot)) 
          stop("either n.tot or n.mc must be specified")
      tii <- estimateGraphFixFast(f.mat, d, q, q.arg, n.mc, n.fast, print.loop.index, 
          ...)
  }
  single.indices <- mainEffectIndices(f.mat=f.mat, d=d, q=q, q.arg=q.arg, n.main=n.main, ...)
  V <- mean(single.indices$V)
  i1 <- as.matrix(round(single.indices$i1, 27))
  rownames(i1) <- paste("X",1:d,sep="")
  tii <- round(tii, 27)
  # estimate cliques
    E <- t(combn(d,2)[,tii[,1] > 0])
    cliques <- maximal.cliques(graph(as.vector(t(E)), d , FALSE))
    tii.scaled <- round(tii[,1,drop=FALSE] / V,27)
    res <- list(d=d, tii=tii, i1=i1, V = V, 
       tii.scaled = tii.scaled, cliques = cliques)
    class(res) <- "graphlist" # define S3 object
  return(res)
}

mainEffectIndices <- function(f.mat, d, q, q.arg, n.main, ...)
{
  X <- matrix(runif(n.main * d), ncol = d)
  for (j in 1:d) X[, j] <- do.call(q[j], c(list(p = X[,j]), q.arg[[j]]))
  Z <- matrix(runif(n.main * d), ncol = d)
  for (j in 1:d) Z[, j] <- do.call(q[j], c(list(p = Z[,j]), q.arg[[j]]))
  yZ <- f.mat(Z,...)
  i1 <- V <- numeric(d)
  names(i1) <- 1:d
  for (i in 1:d){
    ZXi <- Z
    ZXi[,-i] <- X[,-i]
    yZi <- f.mat(ZXi,...)
    i1[i] <- mean(yZi*yZ) - mean(yZi)*mean(yZ) # from Janon et al. (2012, (6))
    V[i] <- mean((yZi^2+yZ^2)/2) - (mean((yZi+yZ)/2))^2
  }
  return(list(i1=i1, V=V))
}