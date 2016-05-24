estShannonf <-
function (X, f) 
{
  X <- as.data.frame(X)
  J <- length(f)
  if (nrow(X) != J) {
    stop("The number of columns in X must be equal to the length of f!")
  }
  if (any(floor(X) != X)) {
    warning("The elements of X should be integers!")
  }
  X <- as.data.frame(X)
  names <- colnames(X)
  ff <- as.factor(f)
  Xs <- split(X, f = ff)
  Xcs <- lapply(X = Xs, FUN = function(x) {
    apply(X = x, MARGIN = 2, FUN = sum)
  })
  tab <- matrix(ncol = ncol(X), nrow = length(Xcs))
  for (i in seq(along.with = Xcs)) {
    tab[i, ] <- Xcs[[i]]
  }
  colnames(tab) <- names
  rownames(tab) <- names(Xcs)
  ngroups <- nrow(tab)
  estimate <- numeric(length = ngroups)
  varest <- numeric(length = ngroups)
  for (i in 1:ngroups) {
    temp <- estShannon(tab[i, ])
    estimate[i] <- temp$estimate
    varest[i] <- temp$varest
  }
  gnames <- rownames(tab)
  names(estimate) <-  names(varest) <- gnames
  out <- list(estimate = estimate, varest = varest)
  return(out)
}

