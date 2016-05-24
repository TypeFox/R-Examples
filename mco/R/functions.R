##
## functions.R - MO test functions
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.uni-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##  Olaf Mersmann    <olafm@statistik.uni-dortmund.de>
##

## nsga2(belegundu, 2, 2, constraints=belegundu.constr, cdim=2, lower.bounds=c(0, 0), upper.bounds=c(5, 3))
belegundu <- function(x) {
  c1 <- 2*x[1]
  c2 <- x[2]
  c(-c1 + c2, c1 + c2)
}

belegundu.constr <- function(x) {
  c(x[1] - x[2] + 1, -x[1] - x[2] + 7)
}

## nsga2(binh1, 2, 2, lower.bounds=c(-5, -5), upper.bounds=c(10, 10))
binh1 <- function(x) {
  c1 <- (x-5)
  c(crossprod(x, x), crossprod(c1, c1))
}

## nsga2(binh2, 2, 2, lower.bounds=c(0, 0), upper.bounds=c(5, 3), constraints=binh2.constr, cdim=2)
binh2 <- function(x) {
  c1 <- 4*x
  c2 <- (x-5)
  c(crossprod(c1, c1), crossprod(c2, c2))
}

binh2.constr <- function(x) {
  c1 <- (x[1] - 5)^2 + x[2]^2 - 25
  c2 <- -(x[1] - 8)^2 - (x[2] + 3)^2 + 7.7
  c(-c1, -c2)
}

## nsga2(binh3, 2, 3, lower.bounds=c(10e-6, 10e-6), upper.bounds=c(10e6, 10e6))
binh3 <- function(x) {
  c1 <- x[1] - 10e6
  c2 <- x[2] - 2e-6
  c3 <- x[1]*x[2] - 2
  c(c1, c2, c3)
}

## FIXME: Fehler in PDF?!
binh4 <- function(x) {
}

## nsga2(deb3, 2, 2, lower.bounds=c(0, 0), upper.bounds=c(1, 1), generations=500)
deb3 <- function(x) {
  g <- function(y) {
    if (y <= 0.4)
      4 - 3*exp(-((y-0.2)/0.02)^2)
    else
      4 - 2*exp(-((y-0.7)/0.2)^2)
  }

  h <- function(a, b) {    
    if (a <= b)
      1 - (a/b)^(0.25 + 3.75*(b - 1))
    else
      0
  }
  ff <- 4*x[1]
  gg <- g(x[2])
  hh <- gg * h(ff, gg)
  c(ff, hh)
}

## nsga2(fonseca1, 2, 2, lower.bounds=c(-100, -100), upper.bounds=c(100, 100))
fonseca1 <- function(x) {
  c1 <- 1 - exp(-(x[1]-1)^2 - (x[2]+1)^2)
  c2 <- 1 - exp(-(x[1]+1)^2 - (x[2]-1)^2)
  c(c1, c2)
}

## nsga2(fonseca2, 2, 2, lower.bounds=c(-4, -4), upper.bounds=c(4, 4))
fonseca2 <- function(x) {
  n <- length(x)
  c1 <- 1/sqrt(n)
  c2 <- x - c1
  c3 <- x + c1
  c4 <- 1 - exp(-crossprod(c2, c2))
  c5 <- 1 - exp(-crossprod(c3, c3))
  c(c4, c5)
}

## nsga2(gianna, 1, 2, lower.bounds=5, upper.bounds=10)
gianna <- function(x) {
  c1 <- 1 / (sqrt(10 - x) + sqrt(x - 5))
  c2 <- 0.04 * (x - 8)^2 + 0.3
  c(c1, c2)
}

## nsga2(hanne1, 2, 2, lower.bounds=c(0, 0), upper.bounds=c(10, 10), constraints=hanne1.constr, cdim=1)
hanne1 <- function(x) { x }

hanne1.constr <- function(x) { sum(x) - 5 }

## nsga2(hanne2, 2, 2, lower.bounds=c(0, 0), upper.bounds=c(10, 10), constraints=hanne2.constr, cdim=1)
hanne2 <- function(x) { x^2 }

hanne2.constr <- function(x) { sum(x) - 5 }

## nsga2(hanne3, 2, 2, lower.bounds=c(0, 0), upper.bounds=c(10, 10), constraints=hanne3.constr, cdim=1)
hanne3 <- function(x) { sqrt(x) }

hanne3.constr <- function(x) { sum(x) - 5 }

## nsga2(hanne4, 2, 2, lower.bounds=c(0, 0), upper.bounds=c(10, 10), constraints=hanne4.constr, cdim=1)
hanne4 <- function(x) { x }

hanne4.constr <- function(x) { x[2] - 5 + 0.5*x[1]*sin(4*x[1]) }

## nsga2(hanne5, 2, 2, lower.bounds=c(0, 0), upper.bounds=c(10, 10), constraints=hanne5.constr, cdim=1)
hanne5 <- function(x) {
  c1 <- 2*pi*(x[2] - trunc(x[2]))
  c2 <- (x[1] - trunc(x[1]))
  c3 <- trunc(x[1]) + 0.5 + c2*sin(c1)
  c4 <- trunc(x[2]) + 0.5 + c2*cos(c1)
  c(c3, c4)
}
hanne5.constr <- function(x) { sum(x) - 5 }

## nsga2(jimenez, 2, 2, lower.bounds=c(0, 0), upper.bounds=c(100, 100), constraints=jimenez.constr, cdim=4)
jimenez <- function(x) { -c(5*x[1] + 3*x[2], 2*x[1] + 8*x[2]) }
jimenez.constr <- function(x) {
  c1 <- x[1] + 4*x[2] - 100
  c2 <- 3*x[1] + 2*x[2] - 150
  c3 <- 200 - 5*x[1] - 3*x[2]
  c4 <- 75 - 2*x[1] - 8*x[2]
  -c(c1, c2, c3, c4)
}

## nsga2(vnt, 2, 3,lower.bounds=rep(-3, 2), upper.bounds=rep(3, 2))
vnt <- function(x) {  
  y <- numeric(3)
  xn <- crossprod(x, x)
  y[1] <- xn/2 + sin(xn);
  y[2] <- (crossprod(c(3, -2), x) + 4)^2/8 + (crossprod(c(1, -1), x) + 1)^2/27 + 15
  y[3] <- 1/(xn + 1) - 1.1*exp(-xn)
  return (y)
}

## nsga2(zdt1, 30, 2, lower.bounds=rep(0, 30), upper.bounds=rep(1, 30))
zdt1 <- function(x) {
  dim <- length(x)
  y1 <- x[1]
  
  g <- 1 + (9 * mean(x[2:dim]))
  y2 <- g * ( 1 - sqrt(y1/g))

  return(c(y1, y2))
}

## nsga2(zdt2, 30, 2, lower.bounds=rep(0, 30), upper.bounds=rep(1, 30))
zdt2 <- function(x) {
  dim <- length(x)
  y1 <- x[1]
  
  g <- 1 + (9 * mean(x[2:dim]))
  y2 <- g * ( 1 - (y1/g)^2)

  return(c(y1, y2))
}

## nsga2(zdt3, 30, 2, lower.bounds=rep(0, 30), upper.bounds=rep(1, 30))
zdt3 <- function(x) {
  dim <- length(x)
  y1 <- x[1]
  
  g <- 1 + (9 * mean(x[2:dim]))
  y2 <- g * ( 1 - sqrt(y1/g) - (y1/g)*sin(10*pi*y1))

  return(c(y1, y2))
}
