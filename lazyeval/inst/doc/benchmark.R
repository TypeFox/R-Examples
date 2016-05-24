## ----, echo = FALSE------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
options(digits = 3)

## ------------------------------------------------------------------------
library(microbenchmark)
library(lazyeval)

microbenchmark(
  substitute = substitute(x),
  lazy = lazy(x)
)

