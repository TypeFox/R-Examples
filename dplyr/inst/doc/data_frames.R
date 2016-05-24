## ---- echo = FALSE, message = FALSE--------------------------------------
knitr::opts_chunk$set(collapse = T, comment = "#>")
options(dplyr.print_min = 4L, dplyr.print_max = 4L)
library(dplyr)

## ------------------------------------------------------------------------
data.frame(x = letters) %>% sapply(class)
data_frame(x = letters) %>% sapply(class)

## ------------------------------------------------------------------------
data_frame(x = 1:3, y = list(1:5, 1:10, 1:20))

## ------------------------------------------------------------------------
data.frame(`crazy name` = 1) %>% names()
data_frame(`crazy name` = 1) %>% names()

## ------------------------------------------------------------------------
data_frame(x = 1:5, y = x ^ 2)

## ------------------------------------------------------------------------
data_frame(x = 1:5) %>% class()

## ------------------------------------------------------------------------
df <- data.frame(a = 1:2, b = 1:2)
str(df[, "a"])

tbldf <- tbl_df(df)
str(tbldf[, "a"])

## ------------------------------------------------------------------------
l2 <- replicate(26, sample(100), simplify = FALSE)
names(l2) <- letters
microbenchmark::microbenchmark(
  as_data_frame(l2),
  as.data.frame(l2)
)

## ------------------------------------------------------------------------
location(iris)

## ------------------------------------------------------------------------
iris2 <- iris
location(iris2)

## ------------------------------------------------------------------------
changes(iris2, iris)

## ------------------------------------------------------------------------
iris2$Sepal.Length <- iris2$Sepal.Length * 2
changes(iris, iris2)

## ------------------------------------------------------------------------
iris3 <- mutate(iris, Sepal.Length = Sepal.Length * 2)
changes(iris3, iris)

