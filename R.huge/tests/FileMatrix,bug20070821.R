library("R.huge")
library("R.utils")

if ("covr" %in% loadedNamespaces())
  options("R.utils::onNonSeekable"="warning")

for (byrow in c(FALSE, TRUE)) {
  pathname <- "example.Rmatrix"
  if (isFile(pathname)) {
    file.remove(pathname)
    if (isFile(pathname))
      stop("File not deleted: ", pathname)
  }

  nrow <- 10
  ncol <- 10
  X <- matrix(0, nrow=nrow, ncol=ncol, byrow=TRUE)
  X[] <- 1:100

  Y <- FileFloatMatrix(pathname, nrow=nrow(X), ncol=ncol(X), byrow=byrow)
  Y[] <- X[]
  stopifnot(identical(X[], Y[]))

  rr <- c(1,3)
  stopifnot(identical(X[rr,], Y[rr,]))

  rr <- 6:2
  cc <- c(1,8,4,2:3)
  X[rr,cc] <- seq(length=length(rr)*length(cc))
  Y[rr,cc] <- X[rr,cc]
  stopifnot(identical(X[], Y[]))

  close(Y)
  rm(Y)
  gc()
}
