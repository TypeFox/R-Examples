## ----, echo = FALSE, message = FALSE-------------------------------------
knitr::opts_chunk$set(comment = "")
library(rjade)

## ------------------------------------------------------------------------
# Compile a Jade template in R
text <- readLines(system.file("examples/test.jade", package = "rjade"))
tpl <- jade_compile(text, pretty = TRUE)

## ------------------------------------------------------------------------
# Render the template
tpl()

## ------------------------------------------------------------------------
tpl(youAreUsingJade = TRUE)

