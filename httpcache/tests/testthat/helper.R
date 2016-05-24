Sys.setlocale("LC_COLLATE", "C") ## What CRAN does
set.seed(999)
options(warn=1)

public <- function (...) with(globalenv(), ...)

public({
    source("helper-mocks.R")
    cacheKeys <- function () ls(envir=httpcache:::cache)
})
