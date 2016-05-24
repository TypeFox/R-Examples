library("globals")

## WORKAROUND: Make sure tests also work with 'covr' package
covr <- ("covr" %in% loadedNamespaces())
if (covr) {
  globalenv <- function() parent.frame()
  baseenv <- function() environment(base::sample)
}

b <- 2
c <- 3
d <- NULL
expr <- substitute({ x <- b; b <- 1; y <- c; z <- d }, env=list())

message("*** findGlobals() ...")

message(" ** findGlobals(..., method='conservative'):")
globalsC <- findGlobals(expr, method="conservative")
print(globalsC)
stopifnot(all(globalsC %in% c("{", "<-", "c", "d")))

message(" ** findGlobals(..., method='liberal'):")
globalsL <- findGlobals(expr, method="liberal")
print(globalsL)
stopifnot(all(globalsL %in% c("{", "<-", "b", "c", "d")))

message(" ** findGlobals(..., method='ordered'):")
globalsI <- findGlobals(expr, method="ordered")
print(globalsI)
stopifnot(all(globalsI %in% c("{", "<-", "b", "c", "d")))

message("*** findGlobals() ... DONE")



message("*** globalsOf() ...")

message(" ** globalsOf(..., method='conservative'):")
globalsC <- globalsOf(expr, method="conservative")
str(globalsC)
stopifnot(all(names(globalsC) %in% c("{", "<-", "c", "d")))
globalsC <- cleanup(globalsC)
str(globalsC)
stopifnot(all(names(globalsC) %in% c("c", "d")))
where <- attr(globalsC, "where")
stopifnot(
  length(where) == length(globalsC),
  identical(where$c, globalenv()),
  identical(where$d, globalenv())
)

message(" ** globalsOf(..., method='liberal'):")
globalsL <- globalsOf(expr, method="liberal")
str(globalsL)
stopifnot(all(names(globalsL) %in% c("{", "<-", "b", "c", "d")))
globalsL <- cleanup(globalsL)
str(globalsL)
stopifnot(all(names(globalsL) %in% c("b", "c", "d")))
where <- attr(globalsL, "where")
stopifnot(
  length(where) == length(globalsL),
  identical(where$b, globalenv()),
  identical(where$c, globalenv()),
  identical(where$d, globalenv())
)

message(" ** globalsOf(..., method='ordered'):")
globalsL <- globalsOf(expr, method="ordered")
str(globalsL)
stopifnot(all(names(globalsL) %in% c("{", "<-", "b", "c", "d")))
globalsL <- cleanup(globalsL)
str(globalsL)
stopifnot(all(names(globalsL) %in% c("b", "c", "d")))
where <- attr(globalsL, "where")
stopifnot(
  length(where) == length(globalsL),
  identical(where$b, globalenv()),
  identical(where$c, globalenv()),
  identical(where$d, globalenv())
)

message("*** globalsOf() ... DONE")


message("*** Subsetting of Globals:")
globalsL <- globalsOf(expr, method="liberal")
globalsS <- globalsL[-1]
stopifnot(length(globalsS) == length(globalsL) - 1L)
stopifnot(identical(class(globalsS), class(globalsL)))
whereL <- attr(globalsL, "where")
whereS <- attr(globalsS, "where")
stopifnot(length(whereS) == length(whereL) - 1L)
stopifnot(identical(whereS, whereL[-1]))


message("*** cleanup() & packagesOf():")
globals <- globalsOf(expr, method="conservative")
str(globals)
stopifnot(all(names(globals) %in% c("{", "<-", "c", "d")))

globals <- as.Globals(globals)
str(globals)
stopifnot(all(names(globals) %in% c("{", "<-", "c", "d")))

globals <- as.Globals(unclass(globals))
str(globals)
stopifnot(all(names(globals) %in% c("{", "<-", "c", "d")))

pkgs <- packagesOf(globals)
print(pkgs)
stopifnot(length(pkgs) == 0L)

globals <- cleanup(globals)
str(globals)
stopifnot(all(names(globals) %in% c("c", "d")))

pkgs <- packagesOf(globals)
print(pkgs)
stopifnot(length(pkgs) == 0L)


message("*** globalsOf() and package functions:")
foo <- globals::Globals
expr <- substitute({ foo(list(a=1)) })
globals <- globalsOf(expr)
str(globals)
stopifnot(all(names(globals) %in% c("{", "foo", "list")))
where <- attr(globals, "where")
stopifnot(
  length(where) == length(globals),
  identical(where$`{`, baseenv()),
  covr || identical(where$foo, globalenv()),
  identical(where$list, baseenv())
)

globals <- cleanup(globals)
str(globals)
stopifnot(all(names(globals) %in% c("foo")))
pkgs <- packagesOf(globals)
stopifnot(pkgs == "globals")


message("*** globalsOf() and core-package functions:")
sample2 <- base::sample
sum2 <- base::sum
expr <- substitute({ x <- sample(10); y <- sum(x); x2 <- sample2(10); y2 <- sum2(x); s <- sessionInfo() }, env=list())
globals <- globalsOf(expr)
str(globals)
stopifnot(all(names(globals) %in% c("{", "<-", "sample", "sample2", "sessionInfo", "sum", "sum2")))
where <- attr(globals, "where")
stopifnot(
  length(where) == length(globals),
  identical(where$`<-`, baseenv()),
  identical(where$sample, baseenv()),
  covr || identical(where$sample2, globalenv())
)

globals <- cleanup(globals)
str(globals)
stopifnot(all(names(globals) %in% c("sample2", "sum2")))
where <- attr(globals, "where")
stopifnot(
  length(where) == length(globals),
  covr || identical(where$sample2, globalenv())
)

globals <- cleanup(globals, drop="primitives")
str(globals)
stopifnot(all(names(globals) %in% c("sample2")))


message("*** globalsOf() - exceptions ...")

rm(list="a")
res <- try({
  globals <- globalsOf({ x <- a }, substitute=TRUE, mustExist=TRUE)
}, silent=TRUE)
stopifnot(inherits(res, "try-error"))

message("*** globalsOf() - exceptions ... DONE")
