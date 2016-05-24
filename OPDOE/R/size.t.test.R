size.t.test <- function(...)
{
  # a ceiling call wrapper to power.t.test
  ceiling(power.t.test(...)$n)
}

delta.t.test <- function(...)
{
  # a ceiling call wrapper to power.t.test
  power.t.test(...)$delta
}
