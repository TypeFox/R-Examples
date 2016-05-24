library("listenv")

ovars <- ls(envir=globalenv())
oopts <- options(warn=1)


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Single-element assignments and subsetting
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
x <- list(a=1, b=2, c=3)
str(x)
y <- as.listenv(x)
print(y)
stopifnot(identical(as.list(y), x))
z <- as.listenv(y)
stopifnot(identical(as.list(y), as.list(z)))

e <- new.env()
e$a <- 1
e$b <- 2
e$c <- 3
y <- as.listenv(e)
print(y)
stopifnot(identical(as.list(y), as.list(e)))

x <- c(a=1, b=2, c=3)
y <- as.listenv(x)
print(y)
stopifnot(identical(as.list(y), as.list(x)))

## Cleanup
options(oopts)
rm(list=setdiff(ls(envir=globalenv()), ovars), envir=globalenv())
