##
## loglognorm.R - Interface to loglognorm.c
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.uni-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##  Olaf Mersmann    <olafm@statistik.uni-dortmund.de>
##

dloglognorm <- function(x, mean=0, sd=1)
  .Call("dloglognorm", x, mean, sd, PACKAGE="loglognorm")

ploglognorm <- function(q, mean=0, sd=1)
  .Call("ploglognorm", q, mean, sd, PACKAGE="loglognorm")

qloglognorm <- function(p, mean=0, sd=1)
  .Call("qloglognorm", p, mean, sd, PACKAGE="loglognorm")

## FIXME: Faster way?
rloglognorm <- function(n, mean=0, sd=1)
  .Call("qloglognorm", runif(n), mean=mean, sd=sd, PACKAGE="loglognorm")

mloglognorm <- function(moment, mean, sd)
  .Call("mloglognorm", moment,  mean, sd, PACKAGE="loglognorm")

eloglognorm <- function(mean, sd)
  .Call("mloglognorm", 1, mean, sd, PACKAGE="loglognorm")

vloglognorm <- function(mean, sd) {
  m1 <- mloglognorm(1, mean, sd)
  m2 <- mloglognorm(2, mean, sd)
  return (m2 - m1^2)
}

