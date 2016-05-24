library("listenv")

ovars <- ls(envir=globalenv())
if (exists("x")) rm(list="x")
if (exists("y")) rm(list="y")

## - - - - - - - - - - - - - - - - - - - - - - - - - -
## Variable in global/parent environment
## - - - - - - - - - - - - - - - - - - - - - - - - - -
message("*** parse_env_subset() on parent environment ...")

target <- parse_env_subset(x, substitute=TRUE)
str(target)
stopifnot(identical(target$envir, environment()),
          target$name == "x", is.na(target$idx), !target$exists)

target <- parse_env_subset("x", substitute=TRUE)
str(target)
stopifnot(identical(target$envir, environment()),
          target$name == "x", is.na(target$idx), !target$exists)

x <- NULL
target <- parse_env_subset(x, substitute=TRUE)
str(target)
stopifnot(identical(target$envir, environment()),
          target$name == "x", is.na(target$idx), target$exists)

target <- parse_env_subset(y, substitute=TRUE)
str(target)
stopifnot(identical(target$envir, environment()),
          target$name == "y", is.na(target$idx), !target$exists)

message("*** parse_env_subset() on parent environment ... DONE")


## - - - - - - - - - - - - - - - - - - - - - - - - - -
## Environment
## - - - - - - - - - - - - - - - - - - - - - - - - - -
message("parse_env_subset() on environment ...")
x <- new.env()

target <- parse_env_subset(x, substitute=TRUE)
str(target)
stopifnot(identical(target$envir, environment()),
          target$name == "x", is.na(target$idx), target$exists)

target <- parse_env_subset(x$a, substitute=TRUE)
str(target)
stopifnot(identical(target$envir, x), target$name == "a", is.na(target$idx), !target$exists)

target <- parse_env_subset("a", envir=x, substitute=TRUE)
str(target)
stopifnot(identical(target$envir, x), target$name == "a", is.na(target$idx), !target$exists)

target <- parse_env_subset(x[["a"]], substitute=TRUE)
str(target)
stopifnot(identical(target$envir, x), target$name == "a", is.na(target$idx), !target$exists)

target <- parse_env_subset("a", envir=x, substitute=TRUE)
str(target)
stopifnot(identical(target$envir, x), target$name == "a", is.na(target$idx), !target$exists)

res <- try(target <- parse_env_subset(1, substitute=FALSE), silent=TRUE)
stopifnot(inherits(res, "try-error"))

res <- try(target <- parse_env_subset(x[[1]], substitute=TRUE), silent=TRUE)
stopifnot(inherits(res, "try-error"))

x$a <- 1
target <- parse_env_subset(x$a, substitute=TRUE)
str(target)
stopifnot(identical(target$envir, x), target$name == "a", is.na(target$idx), target$exists)

message("parse_env_subset() on environment ... DONE")


## - - - - - - - - - - - - - - - - - - - - - - - - - -
## List environment
## - - - - - - - - - - - - - - - - - - - - - - - - - -
message("*** parse_env_subset() on listenv ...")

x <- listenv()

target <- parse_env_subset(x, substitute=TRUE)
str(target)
stopifnot(identical(target$envir, environment()),
          target$name == "x", is.na(target$idx), target$exists)

target <- parse_env_subset(x$a, substitute=TRUE)
str(target)
stopifnot(identical(target$envir, x), target$name == "a", is.na(target$idx), !target$exists)

target <- parse_env_subset(x[["a"]], substitute=TRUE)
str(target)
stopifnot(identical(target$envir, x), target$name == "a", is.na(target$idx), !target$exists)

target <- parse_env_subset("a", envir=x, substitute=TRUE)
str(target)
stopifnot(identical(target$envir, x), target$name == "a", is.na(target$idx), !target$exists)

target <- parse_env_subset(x[[1]], substitute=TRUE)
str(target)
stopifnot(identical(target$envir, x), target$name == "", target$idx == 1, !target$exists)

target <- parse_env_subset(x[[2]], substitute=TRUE)
str(target)
stopifnot(identical(target$envir, x), target$name == "", target$idx == 2, !target$exists)

x$a <- 1
target <- parse_env_subset(x$a, substitute=TRUE)
str(target)
stopifnot(identical(target$envir, x), target$name == "a", target$idx  == 1, target$exists)

target <- parse_env_subset("a", envir=x, substitute=TRUE)
str(target)
stopifnot(identical(target$envir, x), target$name == "a", target$idx  == 1, target$exists)

stopifnot(x$a == 1)
stopifnot(x[[1]] == 1)

target <- parse_env_subset(x[[1]], substitute=TRUE)
str(target)
stopifnot(identical(target$envir, x), target$name == "a", target$idx  == 1, target$exists)


x[[3]] <- 3
target <- parse_env_subset(x[[3]], substitute=TRUE)
str(target)
stopifnot(identical(target$envir, x), target$name == "", target$idx  == 3, target$exists)
stopifnot(x[[3]] == 3)
print(names(x))
stopifnot(identical(names(x), c("a", "", "")))


b <- 1
target <- parse_env_subset(x[[b]], substitute=TRUE)
str(target)
stopifnot(identical(target$envir, x), target$name == "a", target$idx  == 1, target$exists)


x <- listenv()
length(x) <- 2

target <- parse_env_subset(x[[1]], substitute=TRUE)
str(target)
stopifnot(!target$exists)

target <- parse_env_subset(x[[2]], substitute=TRUE)
str(target)
stopifnot(!target$exists)

target <- parse_env_subset(x[[3]], substitute=TRUE)
str(target)
stopifnot(!target$exists)
stopifnot(length(x) == 2)

x[[2]] <- 2
target <- parse_env_subset(x[[2]], substitute=TRUE)
str(target)
stopifnot(target$exists)

x[[4]] <- 4
stopifnot(length(x) == 4)

target <- parse_env_subset(x[[3]], substitute=TRUE)
str(target)
stopifnot(!target$exists)

target <- parse_env_subset(x[1:5], substitute=TRUE)
stopifnot(length(target$idx) == 5, all(target$idx == 1:5))
str(target)

target <- parse_env_subset(x[integer(0L)], substitute=TRUE)
stopifnot(length(target$idx) == 0)
str(target)

target <- parse_env_subset(x[[integer(0L)]], substitute=TRUE)
stopifnot(length(target$idx) == 0)
str(target)

target <- parse_env_subset(x[0], substitute=TRUE)
stopifnot(length(target$idx) == 0)
str(target)

target <- parse_env_subset(x[-1], substitute=TRUE)
stopifnot(length(target$idx) == length(x)-1)
str(target)

## Odds and ends
target <- parse_env_subset(x[[""]], substitute=TRUE)
stopifnot(length(target$idx) == 1L, !target$exists)

message("*** parse_env_subset() on listenv ... DONE")


## - - - - - - - - - - - - - - - - - - - - - - - - - -
## Exception handling
## - - - - - - - - - - - - - - - - - - - - - - - - - -
message("*** parse_env_subset() - exceptions ...")

x <- listenv()

res <- try(target <- parse_env_subset(x[[0]], substitute=TRUE), silent=TRUE)
stopifnot(inherits(res, "try-error"))

res <- try(target <- parse_env_subset("_a", substitute=TRUE), silent=TRUE)
stopifnot(inherits(res, "try-error"))

res <- try(target <- parse_env_subset(1:10, envir=x, substitute=FALSE), silent=TRUE)
stopifnot(inherits(res, "try-error"))

res <- try(target <- parse_env_subset(c("a", "b"), envir=x, substitute=FALSE), silent=TRUE)
stopifnot(inherits(res, "try-error"))

res <- try(target <- parse_env_subset(x@a, substitute=TRUE), silent=TRUE)
stopifnot(inherits(res, "try-error"))

message("*** parse_env_subset() - exceptions ... DONE")


## Cleanup
rm(list=setdiff(ls(envir=globalenv()), ovars), envir=globalenv())
