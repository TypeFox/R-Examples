## ---- echo = FALSE, message = FALSE--------------------------------------
knitr::opts_chunk$set(collapse = T, comment = "#>")
options(dplyr.print_min = 4L, dplyr.print_max = 4L)
library(dplyr)

## ---- eval = FALSE-------------------------------------------------------
#  summarise(per_day, flights = sum(flights))

## ---- eval=FALSE---------------------------------------------------------
#  foo <- function(x) x*x
#  summarise(per_day, flights = foo(sum(flights)) )

