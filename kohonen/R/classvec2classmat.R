classvec2classmat <- function(yvec)
{
  yvec <- factor(yvec)
  nclasses <- nlevels(yvec)

  outmat <- matrix(0, length(yvec), nclasses)
  dimnames(outmat) <- list(NULL, levels(yvec))
  
  for (i in 1:nclasses)
    outmat[which(as.integer(yvec) == i),i] <- 1

  outmat
}

classmat2classvec <- function(ymat, threshold=0)
{
  class.names <- dimnames(ymat)[[2]]
  if (is.null(class.names)) class.names <- 1:ncol(ymat)

  classes <- apply(ymat, 1, function(x) which(x == max(x))[1])
  classes[apply(ymat, 1, max) < threshold] <- NA
  
  class.names[classes]
}
