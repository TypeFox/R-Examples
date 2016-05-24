##
## realistic_df.R - realistic desirability functions
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.tu-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##  Olaf Mersmann    <olafm@statistik.tu-dortmund.de>
##

realisticDF <- function(f, ...)
  UseMethod("realisticDF", f)

## Print method:
print.realistic.desire.function <- function(x, ...) {
  cat("Realistic ")
  print(environment(x)$f)
}
  
realisticDF.default <- function(f, ...)
  stop("Not implemented")

realisticDF.desire.function <- function(f, ...) {  
  ev <- function(x, sd)
    edesire(f, x, sd)
  
  class(ev) <- c("realistic.desire.function", class(f))
  attr(ev, "desire.type") <- paste("Realistic", attr(f, "desire.type"))
  attr(ev, "y.range") <- attr(f, "y.range")
  return(ev) 
}

realisticDF.composite.desire.function <- function(f, ...)
  stop("Please wrap a realistic desirability instead of applying realisticDF to a wrapped desirability.")
