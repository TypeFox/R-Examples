## ----, Setup, echo = FALSE, results = "hide"-----------------------------
set.seed(19790801)
library(assertive)
knitr::opts_chunk$set(error = FALSE)

## ------------------------------------------------------------------------
geomean <- function(x, na.rm = FALSE)
{
  exp(mean(log(x), na.rm = na.rm))
}

## ----, Geomean2----------------------------------------------------------
geomean2 <- function(x, na.rm = FALSE)
{
  assert_is_numeric(x)
  exp(mean(log(x), na.rm = na.rm))
}

## ----, Geomean3----------------------------------------------------------
geomean3 <- function(x, na.rm = FALSE)
{
  assert_is_numeric(x)
  if(!all(is_non_negative(x), na.rm = TRUE)) # Don't worry about NAs here
  {
    warning("x contains negative values, so the geometric mean makes no sense.")
    return(NaN)
  }
  exp(mean(log(x), na.rm = na.rm))
}

## ------------------------------------------------------------------------
geomean4 <- function(x, na.rm = FALSE)
{
  assert_is_numeric(x)
  if(!all(is_non_negative(x), na.rm = TRUE)) # Don't worry about NAs here
  {
    warning("x contains negative values, so the geometric mean makes no sense.")
    return(NaN)
  }
  na.rm <- coerce_to(use_first(na.rm), "logical")
  exp(mean(log(x), na.rm = na.rm))
}

