##
## derringer-suich.R - Derringer-Suich type desirability functions
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.tu-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##  Olaf Mersmann    <olafm@statistik.tu-dortmund.de>
##

##
## The exact distribution for general Derringer-Suich type
## desirabilites is not know. Only special cases can be calculated
## exactly. These are given in Steuer (2005) and summarized in the
## following table:
##
##  Class   | Type                        | Page
## ---------+-----------------------------+------
##  dsLTU11 | (l, t, u, 1, 1)             | 61
##  dsITU11 | (-Inf, t, u, 1, 1)          | c.f. 85
##  dsLTI11 | (l, t, Inf, 1, 1)           | 73
##  dsTTI11 | (t-delta, t, Inf, 1, 1)     | 75
##  dsLTT01 | (l+delta, t, t+delta, 0, 1) | 76
##  dsA1    | (y, d, beta==1)              | 85
##
## Case ITU is a special case of dsA1
##

derringerSuich <- function(y, d, beta) {
  ## Check if 'short' spec is used and convert it into 'long' format.
  if (length(y) == 5 && missing(d) && missing(beta)) {
    d <- c(0, 1, 0)
    beta <- y[4:5]
    y <- y[1:3]
  }

  ## Possibly cast integers to REAL:
  y <- as.numeric(y)
  d <- as.numeric(d)
  beta <- as.numeric(beta)
  
  ev <- function(x, ...)
    .Call("ds_eval", x, y, d, beta)
  
  n <- length(y)
  if (length(d) != n)
    stop("Number of desirabilities does not match number of data points.")
  if (length(beta) != (n-1))
    stop("Number of weights does not match number of data points.")
  if (is.unsorted(y))
    stop("Data points 'y' not ordered.")

  if (any(d < 0 | d > 1))
    stop("Not all desirabilities in the range 0 to 1,")
  if (any(beta <= 0))
    stop("Not all weights are positive.")
  
  class(ev) <- c("derringerSuich", "desire.function")
  ## Check for special cases
  if (length(d) == 3 && all(d == c(0, 1, 0))) { # 'simple' DS
    if (all(beta == 1)) { # (?, ?, ?, 1, 1)
      if (all(is.finite(y))) { # (l, t, u, 1, 1)
        class(ev) <- c("dsLTU11", class(ev))
      } else if (y[3] == Inf) { # (l, t, Inf, 1, 1)
        if (y[1] < y[2]) 
          class(ev) <- c("dsLTI11", class(ev)) # (l, t, Inf, 1, 1)
        else
          class(ev) <- c("dsTTI11", class(ev)) # (t, t, Inf, 1, 1)
      } else if (y[1] == -Inf) { ## Hunch
        class(ev) <- c("dsITU11", class(ev)) # (-Inf, t, u, 1, 1)
      }
    } else if (beta[1] == 0 && beta[2] == 1) { # (?, ?, ?, 0, 1)
      if (y[2] == y[3]) # (l, t, t, 0, 1)
        class(ev) <- c("dsLTT01", class(ev))
    }
  } else if (all(beta == 1)) {
    class(ev) <- c("dsA1", class(ev))
  }

  attr(ev, "desire.type") <- "Derringer-Suich"
  attr(ev, "y.range") <- range(y[is.finite(y)])
  ## Remove unnecessary variables, since they will be saved in ev's environment. 
  rm(n)
  return(ev)
}

## print methods:
print.dsLTU11 <- function(x, ...) {
  e <- environment(x)
  p <- c(e$y[1], e$y[2], e$y[3], 1, 1)
  message("    (", paste(p, collapse=", "),  ") Derringer Suich desirability", sep="")
}

print.dsLTI11 <- function(x, ...) {
  e <- environment(x)
  p <- c(e$y[1], e$y[2], Inf, 1, 1)
  message("    (", paste(p, collapse=", "),  ") Derringer Suich desirability", sep="")
}

print.derringerSuich <- function(x, ...) {
  e <- environment(x)
  message("    Generalized Derringer Suich type desirability")
  message("")
  X <- cbind(format(e$y, width=6),
             format(e$d, width=6),
             c(format(e$beta,width=6), ""))
  colnames(X) <- c("y", "d", "beta")
  print(X, quote=FALSE, right=TRUE)
}

## Case dsLTU11
ddesire.dsLTU11 <- function(x, f, mean=0, sd=1) {
  e <- environment(f)
  .Call("ddsLTU11", x, e$y[1], e$y[2], e$y[3], mean, sd)
}

pdesire.dsLTU11 <- function(q, f, mean=0, sd=1) {
  e <- environment(f)
  .Call("pdsLTU11", q, e$y[1], e$y[2], e$y[3], mean, sd)
}

edesire.dsLTU11 <- function(f, mean=0, sd=1) {
  e <- environment(f)
  .Call("edsLTU11", e$y[1], e$y[2], e$y[3], mean, sd)
}

## Case dsITU11
edesire.dsITU11 <- function(f, mean=0, sd=1) {
  e <- environment(f)
  .Call("edsA1", e$y, c(1, 1, 0), mean, sd);
}
  
## Case dsLTI11
ddesire.dsLTI11 <- function(x, f, mean=0, sd=1) {
  e <- environment(f)
  .Call("ddsLTI11", x, e$y[1], e$y[2], mean, sd)
}

pdesire.dsLTI11 <- function(q, f, mean=0, sd=1) {
  e <- environment(f)
  .Call("pdsLTI11", q, e$y[1], e$y[2], mean, sd)
}

edesire.dsLTI11 <- function(f, mean=0, sd=1) {
  e <- environment(f)
  .Call("edsLTI11", e$y[1], e$y[2], mean, sd)
}

## Case dsA1
## ddesire.dsA1 <- function(x, f, mean=0, sd=1) {
##   e <- environment(f)
##   .Call("ddsA1", x, e$y, e$d, mean, sd);
## }

## pdesire.dsA1 <- function(q, f, mean=0, sd=1) {
##   e <- environment(f)
##   .Call("pdsA1", q, e$y, e$d, mean, sd);
## }

edesire.dsA1 <- function(f, mean=0, sd=1) {
  e <- environment(f)
  .Call("edsA1", e$y, e$d, mean, sd);
}
