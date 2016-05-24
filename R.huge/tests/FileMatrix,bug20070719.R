library("R.huge")
library("R.utils")

if ("covr" %in% loadedNamespaces())
  options("R.utils::onNonSeekable"="warning")

pathname <- "example.Rmatrix"
if (isFile(pathname)) {
  file.remove(pathname)
  if (isFile(pathname))
    stop("File not deleted: ", pathname)
}

nrow <- 50000
ncol <- 20
X <- FileFloatMatrix(pathname, nrow=nrow, ncol=ncol, byrow=TRUE)

n <- 1000
for (ii in 1:(nrow%/%n)) {
  jj <- ((ii-1)*n+1):(ii*n)
  tim <- system.time({
    X[jj,] <- matrix(rnorm(n*ncol), ncol=ncol)
  })[3]
  cat(ii, tim, "s\n")
}

set.seed(1)
for (ii in 1:6) {
  jj <- sample(1:nrow, size=n)
  tim <- system.time({
    y <- X[jj,]
  })[3]
  cat(ii, tim, "s\n")
}

close(X)
delete(X)
