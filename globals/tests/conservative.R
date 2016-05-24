library("globals")

ovars <- ls(envir=globalenv())


## WORKAROUND: Avoid problem reported in testthat Issue #229, which
## causes covr::package_coverage() to given an error. /HB 2015-02-16
suppressWarnings({
  rm(list=c("a", "b", "c", "x", "y", "z", "square",
            "pathname", "url", "filename"))
})


message("Setting up expressions")
exprs <- list(
  A = substitute({ Sys.sleep(1); x <- 0.1 }, env=list()),
  B = substitute({ y <- 0.2 }, env=list()),
  C = substitute({ z <- a+0.3 }, env=list()),
  D = substitute({ pathname <- file.path(dirname(url), filename) }, env=list()),
  E = substitute({ b <- c }, env=list()),
  F = substitute({
    a <- { runif(1) }
    b <- { rnorm(1) }
    x <- a*b; abs(x)
  }, env=list()),
  G = substitute({
    y <- square(a)
  }, env=list()),
  H = substitute({
    b <- a
    a <- 1
  }, env=list())
)

atleast <- list(
  A = c(),
  B = c(),
  C = c("a"),
  D = c("filename"),
  E = c("c"),
  F = c(),
  G = c("a", "square"),
  H = c() ## FIXME: Should be c("a"), cf. Issue #5.
)

not <- list(
  A = c("x"),
  B = c("y"),
  C = c("z"),
  D = c("pathname"),
  E = c("b"),
  F = c("a", "b", "x"),
  G = c(),
  H = c()
)


## Define globals
a <- 3.14
c <- 2.71
square <- function(x) x^2
filename <- "index.html"
# Yes, pretend we forget 'url'

message("Find globals")
for (kk in seq_along(exprs)) {
  key <- names(exprs)[kk]
  expr <- exprs[[key]]
  cat(sprintf("Expression #%d ('%s'):\n", kk, key))
  print(expr)

  names <- findGlobals(expr, method="conservative")
  cat(sprintf("Globals: %s\n", paste(sQuote(names), collapse=", ")))
  stopifnot(all(atleast[[key]] %in% names))
  stopifnot(!any(names %in% not[[key]]))

  globals <- globalsOf(expr, method="conservative")
  cat(sprintf("Globals: %s\n", paste(sQuote(names(globals)), collapse=", ")))
  stopifnot(all(atleast[[key]] %in% names(globals)))
  stopifnot(!any(names(globals) %in% not[[key]]))
  str(globals)

  cat("\n")
}

names <- findGlobals(exprs, method="conservative", unlist=TRUE)
cat(sprintf("Globals: %s\n", paste(sQuote(names), collapse=", ")))


## Cleanup
rm(list=setdiff(ls(envir=globalenv()), ovars), envir=globalenv())
