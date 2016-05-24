library(nlmeU)
dtf <- subset(armd.wide, select = c(visual12, visual24, visual52))
res <-  missPat(dtf, symbols = c("?","+"))

library(testthat)
expect_that(res[1], equals(c("1" = "+??")))