library("listenv")

message("*** undim() ...")

## General
x <- c(a=1, b=2, A=3, B=4)
names <- names(x)
dim <- c(2,2)
dimnames <- list(c("a", "b"), c("A", "B"))

## Basic arrays
y <- array(x, dim=dim, dimnames=dimnames)
names(y) <- names
z <- undim(y)
stopifnot(identical(names(z), names))

## Lists
y <- as.list(x)
dim(y) <- dim
dimnames(y) <- dimnames
names(y) <- names
z <- undim(y)
stopifnot(identical(names(z), names))

## List environments
y <- as.listenv(x)
dim(y) <- dim
dimnames(y) <- dimnames
names(y) <- names
z <- undim(y)
stopifnot(identical(names(z), names))

message("*** undim() ... DONE")
