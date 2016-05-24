estSimpsonf <-
function (X, f)
{
  if (any(floor(X) != X)) {
    warning("The elements of X should be integers!")
  }
  J <- length(f)
  if (nrow(X) != J) {
    stop("The number of columns in X must be equal to the length of f!")
  }
  varestSimpson <- function(x) {
    n <- sum(x)
    estp <- x/n
    S2 <- sum(estp^2)
    S3 <- sum(estp^3)
    2 * (S2 + 2 * (n - 2) * S3 + (3 - 2 * n) * S2^2)/(n *
                                                      (n - 1))
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
  estlist <- lapply(X = Xcs, FUN = function(x) {
    estSimpson(x)
  })
  estv <- unlist(estlist)
  varestlist <- lapply(X = Xcs, FUN = function(x) {
    varestSimpson(x)
  })
  varestv <- unlist(varestlist)
  return(list(estimate = estv, varest = varestv, table = tab))
}

