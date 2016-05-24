##
## harrington2.R - Two sided Harrington type desiraility functions
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.tu-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##  Olaf Mersmann    <olafm@statistik.tu-dortmund.de>
##

harrington2 <- function(LSL, USL, n) {
  ev <- function(y, ...) {
    ys <- (2*y - s)/d
    return(exp(-abs(ys)^n))
  }
  if (USL <= LSL)
    stop("LSL must be smaller than USL.")
  if (n <= 0)
    stop("Invalid exponent 'n'")
  
  ## Precompute:
  s <- USL + LSL
  d <- USL - LSL
  class(ev) <- c("harrington2", "desire.function")
  attr(ev, "desire.type") <- "Two sided Harrington"
  attr(ev, "y.range") <- c(LSL-d, USL + d)
  return(ev)
}

## Print function
print.harrington2 <- function(x, ...) {
  e <- environment(x)
  message("    Two sided Harrington type desirability")
  message("")  
  p <- c(e$LSL, e$USL, e$n)
  names(p) <- c("LSL", "USL", "n")
  message("Parameters:")
  print.default(format(p, width=8), print.gap=2, quote=FALSE, ...)  
}

## Plot function
plot.harrington2 <- function(x, ...) {
  e <- environment(x)
  plot.desire.function(x, ...)
  abline(v=c(e$LSL, e$USL), col="grey", lty=2)
  mtext(c("LSL", "USL"), at=c(e$LSL, e$USL), line=.5, cex=par("cex"))
}

## Density:
dharrington2 <- function(x, LSL, USL, n, mean, sd) {
  s <- USL + LSL
  d <- USL - LSL
  n0 <- 1/n
  n1 <- n0 - 1
  mu.t <- (2/d) * mean - s / d
  sd.t <- (2/d) * sd
  
  c0 <- sqrt(2*pi) * sd.t * x * n
  c1 <- -log(x)
  c1a <- c1^n0
  c1b <- c1^n1
  c2 <- 2*sd.t^2
  c3 <- (c1a - mu.t)^2 / c2
  c4 <- (c1a + mu.t)^2 / c2
  
  r <- c1b/c0 * (exp(-c3) + exp(-c4))
  return(r)  
}

ddesire.harrington2 <- function(x, f, mean, sd) {
  e <- environment(f)
  dharrington2(x, e$LSL, e$USL, e$n, mean, sd)
}

## CDF
pharrington2 <- function(q, LSL, USL, n, mean, sd) {
  s <- USL + LSL
  d <- USL - LSL

  mu.t <- 2 / d * mean - s / d
  sd.t <- (2 / d) * sd

  n0 <- 1/n
  c1 <- -log(q)^n0
            
  return(2 - pnorm((c1 - mu.t)/sd.t) - pnorm((c1 + mu.t)/sd.t))
}

pdesire.harrington2 <- function(q, f, mean, sd) {
  e <- environment(f)
  pharrington2(q, e$LSL, e$USL, e$n, mean, sd)
}

## Quantiles
qharrington2 <- function(p, LSL, USL, n, mean, sd) {
  f <- function(q) 
    (p - pharrington2(q, LSL, USL, n, mean, sd))^2
  ## FIXME: Ugly search
  optimize(f, c(0, 1))$minimum
}

qdesire.harrington2 <- function(p, f, mean, sd) {
  e <- environment(f)
  qharrington2(p, e$LSL, e$USL, e$n, mean, sd)
}

## Random numbers
rharrington2 <- function(ns, LSL, USL, n, mean, sd)
  harrington2(LSL, USL, n)(rnorm(ns, mean, sd))

## Expectation
eharrington2 <- function(LSL, USL, n, mean, sd) 
  stop("Not implemented")

## Variance
vharrington2 <- function(LSL, USL, n, mean, sd) 
  stop("Not implemented")

