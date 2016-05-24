library("listenv")

ovars <- ls(envir=globalenv())
oopts <- options(warn=1)
map <- listenv:::map

message("* get_variable() - multi-dimensional list environments ...")

x <- listenv()
length(x) <- 6
dim(x) <- c(2,3)

for (ii in seq_along(x)) {
  stopifnot(is.null(x[[ii]]))
  idx <- arrayInd(ii, .dim=dim(x))
  stopifnot(is.null(x[[idx[1],idx[2]]]))
  varV <- get_variable(x, ii, create=FALSE)
  varA <- get_variable(x, idx, create=FALSE)
  stopifnot(identical(varA, varV))
}

x[1:6] <- 1:6
for (ii in seq_along(x)) {
  stopifnot(identical(x[[ii]], ii))
  idx <- arrayInd(ii, .dim=dim(x))
  stopifnot(identical(x[[idx[1],idx[2]]], ii))

  varV <- get_variable(x, ii)
  varA <- get_variable(x, idx)
  stopifnot(identical(varA, varV))
}


message("* get_variable() - multi-dimensional list environments ... DONE")

## Cleanup
options(oopts)
rm(list=setdiff(ls(envir=globalenv()), ovars), envir=globalenv())
