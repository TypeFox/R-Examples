## ------------------------------------------------------------------------
fetch_hook_gh_description <- function(key, namespace) {
  if (!isTRUE(unname(capabilities("libcurl")))) {
    stop("This vignette requires libcurl support in R to run")
  }
  fmt <- "https://raw.githubusercontent.com/%s/master/DESCRIPTION"
  path <- tempfile("gh_description_")
  on.exit(file.remove(path))
  code <- download.file(sprintf(fmt, key), path, method="libcurl")
  if (code != 0L) {
    stop("Error downloading file")
  }
  as.list(read.dcf(path)[1, ])
}

## ------------------------------------------------------------------------
st <- storr::storr_external(storr::driver_environment(),
                            fetch_hook_gh_description)

## ------------------------------------------------------------------------
st$list()

## ------------------------------------------------------------------------
d <- st$get("richfitz/storr")

## ------------------------------------------------------------------------
identical(st$get("richfitz/storr"), d)

## ------------------------------------------------------------------------
st$list()

## ------------------------------------------------------------------------
tryCatch(st$get("richfitz/no_such_repo"),
         KeyErrorExternal=function(e)
           message(sprintf("** Repository %s not found", e$key)))

## ------------------------------------------------------------------------
st_rds <- st$export(storr::storr_rds(tempfile(), mangle_key=TRUE))
st_rds$list()
st_rds$get("richfitz/storr")$Version

## ------------------------------------------------------------------------
st_rds$destroy()

## ------------------------------------------------------------------------
f <- function(a, b) {
  message(sprintf("Computing f(%.3f, %.3f)", a, b))
  ## ...expensive computation here...
  list(a, b)
}

## ------------------------------------------------------------------------
pars <- data.frame(id=as.character(1:10), a=runif(10), b=runif(10),
                   stringsAsFactors=FALSE)

## ------------------------------------------------------------------------
hook <- function(key, namespace) {
  p <- pars[match(key, pars$id), -1]
  f(p$a, p$b)
}

st <- storr::storr_external(storr::driver_environment(), hook)

## ------------------------------------------------------------------------
x <- st$get("1")

## ------------------------------------------------------------------------
identical(st$get("1"), x)

## ------------------------------------------------------------------------
st <- storr::storr_environment()
st$set("experiment1", pars, namespace="parameters")
st$set("experiment1", f, namespace="functions")

hook2 <- function(key, namespace) {
  f <- st$get(namespace, namespace="functions")
  pars <- st$get(namespace, namespace="parameters")
  p <- pars[match(key, pars$id), -1]
  f(p$a, p$b)
}

st_use <- storr::storr_external(st$driver, hook2)

x1 <- st_use$get("1", "experiment1")
x2 <- st_use$get("1", "experiment1")

## ------------------------------------------------------------------------
memoise <- function(f, driver=storr::driver_environment()) {
  force(f)
  st <- storr::storr(driver)
  function(...) {
    key <- digest::digest(list(...))
    tryCatch(
      st$get(key),
      KeyError=function(e) {
        ans <- f(...)
        st$set(key, ans)
        ans
      })
  }
}

## ------------------------------------------------------------------------
f <- function(x) {
  message("computing...")
  x * 2
}

## ------------------------------------------------------------------------
g <- memoise(f)

## ------------------------------------------------------------------------
g(1)

## ------------------------------------------------------------------------
g(1)

## ----eval=FALSE, echo=FALSE----------------------------------------------
#  f <- function(x) x * 2
#  g <- memoise(f)
#  h <- memoise::memoise(f)
#  microbenchmark::microbenchmark(f(1), g(1), h(1))

