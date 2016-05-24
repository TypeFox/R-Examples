##
## truncnorm.R - Interface to truncnorm.c
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.uni-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##  Olaf Mersmann    <olafm@statistik.uni-dortmund.de>
##

dtruncnorm <- function(x, a=-Inf, b=Inf, mean=0, sd=1)
  .Call("do_dtruncnorm", x, a, b, mean, sd)

ptruncnorm <- function(q, a=-Inf, b=Inf, mean=0, sd=1)
  .Call("do_ptruncnorm", q, a, b, mean, sd)

qtruncnorm <- function(p, a=-Inf, b=Inf, mean=0, sd=1)
  .Call("do_qtruncnorm", p, a, b, mean, sd)

rtruncnorm <- function(n, a=-Inf, b=Inf, mean=0, sd=1) {
  if (length(n) > 1)
    n <- length(n)
  if (length(n) > 1)
    n <- length(n)
  else if (!is.numeric(n))
    stop("non-numeric argument n.")
  .Call("do_rtruncnorm", as.integer(n), a, b, mean, sd)
}

etruncnorm <- function(a=-Inf, b=Inf, mean=0, sd=1)
  .Call("do_etruncnorm", a, b, mean, sd)

vtruncnorm <- function(a=-Inf, b=Inf, mean=0, sd=1)
  .Call("do_vtruncnorm", a, b, mean, sd)
