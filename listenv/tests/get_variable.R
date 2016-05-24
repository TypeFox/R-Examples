library("listenv")

ovars <- ls(envir=globalenv())
oopts <- options(warn=1)
map <- listenv:::map

x <- listenv()
length(x) <- 3L
names(x) <- c("a", "b", "c")
stopifnot(length(x) == 3L)
print(map(x))

var <- get_variable(x, "a")
stopifnot(!is.na(var))
stopifnot(length(x) == 3L)
print(map(x))

var <- get_variable(x, "b")
stopifnot(!is.na(var))
stopifnot(length(x) == 3L)
print(map(x))

var <- get_variable(x, "c")
stopifnot(!is.na(var))
stopifnot(length(x) == 3L)
print(map(x))

var <- get_variable(x, "d")
stopifnot(!is.na(var))
stopifnot(length(x) == 4L)
print(map(x))

var <- get_variable(x, 4L)
stopifnot(!is.na(var))
stopifnot(length(x) == 4L)
print(map(x))

x$b <- 2
var <- get_variable(x, "b")
stopifnot(!is.na(var))
stopifnot(length(x) == 4L)
print(map(x))

var <- get_variable(x, length(x) + 1L)
stopifnot(length(x) == 5L)
print(names(x))
print(map(x))

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Allocation
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
x <- listenv()
length(x) <- 3L
print(x[[1]])
print(x[[2]])
print(x[[3]])

## Out-of-bound subsetting
res <- try(x[[0]], silent=TRUE)
stopifnot(inherits(res, "try-error"))

## Out-of-bound subsetting
res <- try(x[[4]], silent=TRUE)
stopifnot(inherits(res, "try-error"))

print(get_variable(x, 1L, mustExist=FALSE))
print(get_variable(x, 2L, mustExist=FALSE))
print(get_variable(x, 3L, mustExist=FALSE))

## Out-of-bound element
res <- try(var <- get_variable(x, 0L, mustExist=TRUE), silent=TRUE)
stopifnot(inherits(res, "try-error"))

## Out-of-bound element
res <- try(var <- get_variable(x, length(x) + 1L, mustExist=TRUE), silent=TRUE)
stopifnot(inherits(res, "try-error"))


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Exception handling
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
x <- listenv()
length(x) <- 3L
names(x) <- c("a", "b", "c")

## Non-existing element
res <- try(var <- get_variable(x, "z", mustExist=TRUE), silent=TRUE)
stopifnot(inherits(res, "try-error"))

res <- try(var <- get_variable(x, c("a", "b")), silent=TRUE)
stopifnot(inherits(res, "try-error"))

res <- try(var <- get_variable(x, 1+2i), silent=TRUE)
stopifnot(inherits(res, "try-error"))



## Cleanup
options(oopts)
rm(list=setdiff(ls(envir=globalenv()), ovars), envir=globalenv())
