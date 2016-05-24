##
## dfun.R - general desirability function utilities
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.tu-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##  Olaf Mersmann    <olafm@statistik.tu-dortmund.de>
##

##
## Desirability function objects:
##
## All desirabilities are native R functions. They may have several
## attributes. These currently include:
##
## * y.range - Range to be plotted
## * desire.type - Pretty name for printing
##

## Useful utility:
is.desirability <- function(x) 
  is.function(x) & "desire.function" %in% class(x)

is.composite.desirability <- function(x)
  is.function(x) & "composite.desire.function" %in% class(x)

## Default print method
print.desire.function <- function(x, ...) 
  message(attr(x, "desire.type"), " type desirability")

## Default plot method
plot.desire.function <- function(x, n=600,
                                 xlim=NULL, ylim=c(0, 1),
                                 xlab="Value", ylab="Desirability",
                                 ..., main) {
  if (is.null(xlim))
    xlim <- attr(x, "y.range")
  plot.new()
  plot.window(xlim, ylim)
  box(); axis(1); axis(2)
  if (missing(main))
    main <- paste("Type: ", attr(x, "desire.type"))
  title(main=main, xlab=xlab, ylab=ylab)

  ## Don't use xlim as range. Use 'real' range of x axis.
  xrng <- par("usr")[1:2]
  z <- seq(xrng[1], xrng[2], length.out=n)
  y <- x(z)
  abline(h=c(0, 1), col="grey", lty=2)
  lines(z, y, ...)
}

plot.realistic.desire.function <- function(x, sd, n=600, ...,
                                           xlim=NULL, ylim=c(0, 1),
                                           xlab="Value", ylab="Desirability") {
  if (is.null(xlim))
    xlim <- attr(x, "y.range")
  plot.new()
  plot.window(xlim, ylim)
  box(); axis(1); axis(2)
  title(main=paste("Type: ", attr(x, "desire.type")),
        xlab=xlab, ylab=ylab)

  ## Don't use xlim as range. Use 'real' range of x axis.
  xrng <- par("usr")[1:2]
  z <- seq(xrng[1], xrng[2], length.out=n)
  y <- x(z, sd)
  abline(h=c(0, 1), col="grey", lty=2)
  lines(z, y, ...)  
}

## Distribution related functions:
ddesire <- function(x, f, mean=0, sd=1)
  UseMethod("ddesire", f)

pdesire <- function(q, f, mean=0, sd=1)
  UseMethod("pdesire", f)

qdesire <- function(p, f, mean=0, sd=1)
  UseMethod("qdesire", f)

rdesire <- function(n, f, mean=0, sd=1)
  UseMethod("rdesire", f)

edesire <- function(f, mean=0, sd=1)
  UseMethod("edesire", f)

vdesire <- function(f, mean=0, sd=1)
  UseMethod("vdesire", f)

##
## These functions try to provide sane default behavior in case we
## have no extra knowledge about the desirabilities distribution.
## This revolves around the generation of a sample from the
## desirabilities distribution and then estimating the desired
## property
##

ddesire.default <- function(x, f, mean=0, sd=1) 
  stop("Not implemented.")

pdesire.default <- function(q, f, mean=0, sd=1) {
  warning("Using finite sample estimation.")
  s <- rdesire(100000, f, mean, sd)
  r <- sapply(q, function(qq) sum(s < qq)/100000)
  return(r)
}

qdesire.default <- function(p, f, mean=0, sd=1) {
  warning("Using finite sample estimation.")
  s <- rdesire(100000, f, mean, sd)
  return(quantile(s, p))
}

rdesire.default <- function(n, f, mean=0, sd=1) {
  return(f(rnorm(n, mean, sd)))
}

edesire.default <- function(f, mean=0, sd=1) { 
  warning("Using finite sample estimation.")
  n <- max(length(mean), length(sd))
  mean <- rep(mean, length.out=n)
  sd <- rep(sd, length.out=n)
  res <- sapply(1:n, function(i) mean(rdesire(100000, f, mean[i], sd[i])))
  return(res)
}

vdesire.default <- function(f, mean=0, sd=1) { 
  warning("Using finite sample estimation.")  
  s <- rdesire(100000, f, mean, sd)  
  return(var(s))
}

