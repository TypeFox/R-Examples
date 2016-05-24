##
## harrington1.R - One sided  Harrington type desiraility functions
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.tu-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##  Olaf Mersmann    <olafm@statistik.tu-dortmund.de>
##

h1.solve.params <- function(y1, d1, y2, d2) {
  ## Solve for constants b0, b1
  ## See Trautmann Diss. p. 13-14
  X <- cbind(1, c(y1, y2))
  b <- solve(X, -log(-log(c(d1, d2))))
  return(b)
}

harrington1 <- function(y1, d1, y2, d2) {
  ev <- function(y, ...) {
    ys <- b0 + b1*y
    return(exp(-exp(-ys)))
  }
  b <- h1.solve.params(y1, d1, y2, d2)
  b0 <- b[1];  b1 <- b[2]
  
  class(ev) <- c("harrington1", "desire.function")
  attr(ev, "desire.type") <- "One sided Harrington"
  attr(ev, "y.range") <- c(y1, y2)
  ## Remove cruft to save space
  rm(b)
  return(ev)
}

## print method
print.harrington1 <- function(x, ...) {
  e <- environment(x)
  message("    One sided Harrington type desirability")
  message("")  
  pi <- c(e$y1, e$d1, e$y2, e$d2)
  names(pi) <- c("y1", "d1", "y2", "d2")
  pc <- c(e$b0, e$b1)
  names(pc) <- c("b0", "b1")
  message("Input parameters:")
  print.default(format(pi, width=8), print.gap=2, quote=FALSE, ...)
  message("Computed parameters:")
  print.default(format(pc, width=8), print.gap=2, quote=FALSE, ...)  
}

## Density
dharrington1 <- function(x, y1, d1, y2, d2, mean, sd) {
  b <- h1.solve.params(y1, d1, y2, d2)
  mu.t <- -(b[1] + b[2] * mean)
  ## OME: sigma.t^2 = b[2]^2 * sd^2 => abs!
  sd.t <- abs(b[2]*sd)
  dloglognorm(x, mu.t, sd.t)
}

ddesire.harrington1 <- function(x, f, mean, sd) {
  e <- environment(f)
  mu.t <- -(e$b0 + e$b1*mean)
  ## OME: sigma.t^2 = e$b1^2 * sd^2 => abs!
  sd.t <- abs(e$b1*sd)
  dloglognorm(x, mu.t, sd.t)
}

## CDF
pharrington1 <- function(q, y1, d1, y2, d2, mean, sd) {
  b <- h1.solve.params(y1, d1, y2, d2)
  mu.t <- -(b[1] + b[2] * mean)
  sd.t <- b[2]*sd
  ploglognorm(q, mu.t, sd.t)
}

pdesire.harrington1 <- function(q, f, mean, sd) {
  e <- environment(f)
  mu.t <- -(e$b0 + e$b1*mean)
  sd.t <- e$b1*sd
  ploglognorm(q, mu.t, sd.t)
}

## Quantiles
qharrington1 <- function(p, y1, d1, y2, d2, mean, sd) {
  b <- h1.solve.params(y1, d1, y2, d2)
  mu.t <- -(b[1] + b[2] * mean)
  sd.t <- b[2]*sd
  qloglognorm(p, mu.t, sd.t)
}

qdesire.harrington1 <- function(p, f, mean, sd) {
  e <- environment(f)
  mu.t <- -(e$b0 + e$b1*mean)
  sd.t <- e$b1*sd
  qloglognorm(p, mu.t, sd.t)
}

## Random numbers
rharrington1 <- function(n, y1, d1, y2, d2, mean, sd)
  harrington1(y1, d1, y2, d2)(rnorm(n, mean, sd))

## Expectation
eharrington1 <- function(y1, d1, y2, d2, mean, sd) {
  b <- h1.solve.params(y1, d1, y2, d2)
  mu.t <- -(b[1] + b[2] * mean)
  sd.t <- b[2]*sd
  eloglognorm(mu.t, sd.t)
}

edesire.harrington1 <- function(f, mean, sd) {
  e <- environment(f)
  mu.t <- -(e$b0 + e$b1*mean)
  sd.t <- e$b1*sd
  eloglognorm(mu.t, sd.t)
}

## Variance
vharrington1 <- function(y1, d1, y2, d2, mean, sd) {
  b <- h1.solve.params(y1, d1, y2, d2)
  mu.t <- -(b[1] + b[2] * mean)
  sd.t <- b[2]*sd
  vloglognorm(mu.t, sd.t)
}

vdesire.harrington1 <- function(f, mean, sd) {
  e <- environment(f)
  mu.t <- -(e$b0 + e$b1*mean)
  sd.t <- e$b1*sd
  vloglognorm(mu.t, sd.t)
}
