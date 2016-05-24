## ---- echo = FALSE, message = FALSE--------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
options(tibble.print_min = 4L, tibble.print_max = 4L)
library(tibble)

## ------------------------------------------------------------------------
data_frame(x = letters)

## ------------------------------------------------------------------------
data_frame(x = 1:3, y = list(1:5, 1:10, 1:20))

## ------------------------------------------------------------------------
names(data.frame(`crazy name` = 1))
names(data_frame(`crazy name` = 1))

## ------------------------------------------------------------------------
data_frame(x = 1:5, y = x ^ 2)

## ------------------------------------------------------------------------
l <- replicate(26, sample(100), simplify = FALSE)
names(l) <- letters

microbenchmark::microbenchmark(
  as_data_frame(l),
  as.data.frame(l)
)

## ------------------------------------------------------------------------
data_frame(x = 1:1000)

## ------------------------------------------------------------------------
df1 <- data.frame(x = 1:3, y = 3:1)
class(df1[, 1:2])
class(df1[, 1])

df2 <- data_frame(x = 1:3, y = 3:1)
class(df2[, 1:2])
class(df2[, 1])

## ------------------------------------------------------------------------
class(df2[[1]])
class(df2$x)

## ---- error = TRUE-------------------------------------------------------
df <- data.frame(abc = 1)
df$a

df2 <- data_frame(abc = 1)
df2$a

