## ----, Setup, echo = FALSE, results = "hide"-----------------------------
set.seed(19790801)
library(assertive)
knitr::opts_chunk$set(error = FALSE)

## ----, IsNumeric---------------------------------------------------------
is_numeric(1:6)
is_numeric(letters)

## ----, IsNonNegative-----------------------------------------------------
is_non_negative(rnorm(6))

