### R code from vignette source 'bvpTests.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: preliminaries
###################################################
library("bvpSolve")
options(prompt = " ")
options(continue = "  ")
options(width=70)


###################################################
### code chunk number 2: bvpTests.Rnw:165-168
###################################################
Prob1 <- function(t, y, pars) {
   list(c( y[2] , y[1]/xi ))
}


###################################################
### code chunk number 3: bvpTests.Rnw:171-181
###################################################
xi <- 0.1
print(system.time(
  mod1 <- bvpshoot(yini = c(1, NA), yend = c(0, NA), x = seq(0, 1, by=0.01),
            func = Prob1, guess = 0)))
print(system.time(
  mod1 <- bvptwp(yini = c(1, NA),yend=c(0, NA),x = seq(0, 1, by = 0.01),
            func = Prob1)))
print(system.time(
  mod1 <- bvpcol(yini = c(1, NA),yend=c(0, NA),x = seq(0, 1, by = 0.01),
            func = Prob1)))


###################################################
### code chunk number 4: bvpTests.Rnw:184-187
###################################################
xi <-0.01
mod2  <- bvptwp(yini = c(1, NA), yend = c(0, NA), x = seq(0, 1, by = 0.01),
            func = Prob1)


###################################################
### code chunk number 5: bvpTests.Rnw:190-193
###################################################
xi <-0.001
mod3  <- bvptwp(yini = c(1, NA), yend = c(0, NA), x = seq(0, 1, by = 0.01),
           func = Prob1)


###################################################
### code chunk number 6: prob1
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 1")

# exact solution
curve(exp(-x/sqrt(xi))-exp((x-2)/sqrt(xi))/(1-exp(-2/sqrt(xi))),
      0, 1, add = TRUE, type = "p")


###################################################
### code chunk number 7: prob1
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 1")

# exact solution
curve(exp(-x/sqrt(xi))-exp((x-2)/sqrt(xi))/(1-exp(-2/sqrt(xi))),
      0, 1, add = TRUE, type = "p")


###################################################
### code chunk number 8: bvpTests.Rnw:221-227
###################################################
Prob2 <- function(t, y, pars) {
  list( y[2]/xi )
}
xi <-0.2
mod1 <- bvpshoot(yini = c(1, NA), yend = c(0, NA), x = seq(0, 1, by = 0.01),
                  order = 2, func = Prob2, guess = 0)


###################################################
### code chunk number 9: bvpTests.Rnw:231-240
###################################################
xi <-0.1
mod2 <- bvptwp(yini = c(1, NA), yend = c(0, NA), x = seq(0, 1, by = 0.01),
                order = 2, func = Prob2, atol = 1e-10)
xi <- 0.01
mod3 <- bvptwp(yini = c(1, NA), yend = c(0, NA), x = seq(0, 1, by = 0.01),
                order = 2, func = Prob2, atol = 1e-10)
xi <- 0.001
mod4 <- bvpcol(yini = c(1, NA), yend = c(0, NA), x = seq(0, 1, by = 0.01),
                order = 2, func = Prob2, atol = 1e-10)


###################################################
### code chunk number 10: prob2
###################################################
plot(mod1, mod2, mod3, mod4, which = 1, lty = 1, main = "test problem 2")
xi <- 0.01
curve((1-exp((x-1)/xi))/(1-exp(-1/xi)), 0, 1, type = "p", add = TRUE)


###################################################
### code chunk number 11: prob2
###################################################
plot(mod1, mod2, mod3, mod4, which = 1, lty = 1, main = "test problem 2")
xi <- 0.01
curve((1-exp((x-1)/xi))/(1-exp(-1/xi)), 0, 1, type = "p", add = TRUE)


###################################################
### code chunk number 12: bvpTests.Rnw:266-273
###################################################
Prob3 <- function(x, y, pars) {
  list(c( y[2],
         1/xi * (-(2+cos(pi*x)) * y[2] + y[1]-
          (1 + xi*pi*pi) * cos(pi*x)-
          (2 + cos(pi*x))* pi * sin(pi*x))
      ))
}


###################################################
### code chunk number 13: bvpTests.Rnw:275-278
###################################################
xi <-0.1
mod1 <- bvpshoot(yini = c(-1, NA), yend = c(-1, NA), 
                   x = seq(-1, 1, by=0.01), func = Prob3, guess = 0)


###################################################
### code chunk number 14: bvpTests.Rnw:280-286
###################################################
xi <-0.01
mod2 <- bvptwp(yini = c(-1, NA), yend = c(-1, NA), 
                 x = seq(-1, 1, by=0.01), func = Prob3)
xi <-0.001
mod3 <- bvpcol(yini = c(-1, NA), yend = c(-1, NA), 
                 x = seq(-1, 1, by=0.01), func = Prob3)


###################################################
### code chunk number 15: prob3
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 3")
curve(cos(pi*x), -1, 1, type = "p", add = TRUE)


###################################################
### code chunk number 16: prob3
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 3")
curve(cos(pi*x), -1, 1, type = "p", add = TRUE)


###################################################
### code chunk number 17: bvpTests.Rnw:311-315
###################################################
Prob4 <- function(t, y, pars) {
  list((-y[2] + (1+xi)*y[1])/xi )
}
yini <- c(1 + exp(-2), NA)


###################################################
### code chunk number 18: bvpTests.Rnw:317-329
###################################################
xi   <- 0.5
yend <- c(1 + exp(-2*(1+xi)/xi), NA)
mod1  <- bvpshoot(yini = yini, yend = yend, x = seq(-1, 1, by = 0.01),
                   order = 2, func = Prob4, guess = 0)
xi   <- 0.1
yend <- c(1 + exp(-2*(1+xi)/xi), NA)
mod2  <- bvptwp(yini = yini, yend = yend, x = seq(-1, 1, by = 0.01),
                 order = 2, func = Prob4)
xi <- 0.01
yend <- c(1 + exp(-2*(1+xi)/xi), NA)
mod3  <- bvptwp(yini = yini,yend = yend,x = seq(-1, 1, by = 0.01),
                  order = 2, func = Prob4)


###################################################
### code chunk number 19: prob4
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 4")
curve(exp(x-1) + exp(-(1+xi)*(1+x)/xi), -1, 1, type = "p", add = TRUE)


###################################################
### code chunk number 20: prob4
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 4")
curve(exp(x-1) + exp(-(1+xi)*(1+x)/xi), -1, 1, type = "p", add = TRUE)


###################################################
### code chunk number 21: bvpTests.Rnw:354-358
###################################################
Prob5 <- function(x, y, pars) {
  list(c( y[2],
          x * y[2] + y[1] - (1+pi*pi) * cos(pi*x) + pi*x*sin(pi*x) ))
}


###################################################
### code chunk number 22: bvpTests.Rnw:360-363
###################################################
xi <- 0.1
mod1 <- bvpshoot(yini = c(-1, NA), yend = c(-1, NA),
                   x=seq(-1, 1, by = 0.01), func = Prob5, guess = 0)


###################################################
### code chunk number 23: bvpTests.Rnw:365-367
###################################################
mod2 <- bvptwp(yini = c(-1, NA), yend = c(-1, NA),
                 x=seq(-1, 1, by = 0.01), func = Prob5)


###################################################
### code chunk number 24: prob5
###################################################
plot(mod1, mod2, which = 1, lty = 1, main = "test problem 5")
curve(cos(pi*x), -1, 1, type = "p", add = TRUE)


###################################################
### code chunk number 25: prob5
###################################################
plot(mod1, mod2, which = 1, lty = 1, main = "test problem 5")
curve(cos(pi*x), -1, 1, type = "p", add = TRUE)


###################################################
### code chunk number 26: bvpTests.Rnw:396-399
###################################################
Prob6 <- function(t, y, pars) {
  list(1/xi * (-t*y[2] - xi*pi*pi*cos(pi*t) - pi*t*sin(pi*t)) )
}


###################################################
### code chunk number 27: bvpTests.Rnw:401-410
###################################################
xi    <- 0.1
mod1 <- bvpshoot(yini = c(-2, NA), yend = c(0, NA), 
            order = 2, x = seq(-1, 1, by = 0.01), func = Prob6, guess = 0)
xi    <- 0.01
mod2 <- bvptwp(yini = c(-2, NA), yend = c(0, NA), 
            order = 2, x = seq(-1, 1, by = 0.01), func = Prob6)
xi    <- 0.001
mod3 <- bvptwp(yini = c(-2, NA), yend = c(0, NA), 
            order = 2, x = seq(-1, 1, by = 0.01), func = Prob6)


###################################################
### code chunk number 28: prob6
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 6")
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
curve(cos(pi*x) + erf(x/sqrt(2*xi))/erf(1/sqrt(2*xi)), -1, 1, 
      type = "p", add = TRUE)



###################################################
### code chunk number 29: prob6
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 6")
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
curve(cos(pi*x) + erf(x/sqrt(2*xi))/erf(1/sqrt(2*xi)), -1, 1, 
      type = "p", add = TRUE)



###################################################
### code chunk number 30: bvpTests.Rnw:438-444
###################################################
prob7 <- function(x, y, pars) {
  list(c(y[2],
     1/xi * (-x*y[2]+y[1] - (1+xi*pi*pi)*cos(pi*x)-pi*x*sin(pi*x)))   )
}

x  <- seq(-1, 1, by = 0.01)


###################################################
### code chunk number 31: bvpTests.Rnw:446-450
###################################################
xi   <- 0.01
mod1  <- bvptwp(yini = c(-1, NA), yend = c(1, NA), x = x, func = prob7)
xi   <- 0.001
mod2 <- bvptwp(yini = c(-1, NA), yend = c(1, NA), x = x, func = prob7)


###################################################
### code chunk number 32: bvpTests.Rnw:453-456
###################################################
xi <- 0.0005
mod3  <- bvptwp(yini = c(-1, NA), yend = c(1, NA), x = x, func = prob7,
  xguess = mod2[,1], yguess = t(mod2[,-1]))


###################################################
### code chunk number 33: prob7
###################################################
plot(mod1, mod2, mod3, which = c(2,1), type = "l", lty = 1,
     main = c("dy", "y"), xlab = "x", ylab = "y")
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
curve(cos(pi*x) + x + (x*erf(x/sqrt(2*xi))+sqrt(2*xi/pi)*exp(-x^2/2/xi))/
         (erf(1/(2*xi))+sqrt(2*xi/pi)*exp(-1/2/xi)),
         -1, 1, type = "p", add = TRUE)


###################################################
### code chunk number 34: prob7
###################################################
plot(mod1, mod2, mod3, which = c(2,1), type = "l", lty = 1,
     main = c("dy", "y"), xlab = "x", ylab = "y")
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
curve(cos(pi*x) + x + (x*erf(x/sqrt(2*xi))+sqrt(2*xi/pi)*exp(-x^2/2/xi))/
         (erf(1/(2*xi))+sqrt(2*xi/pi)*exp(-1/2/xi)),
         -1, 1, type = "p", add = TRUE)


###################################################
### code chunk number 35: bvpTests.Rnw:486-491
###################################################
prob8 <- function(x, y, pars) {
  list(-1/xi*y[2])
}

x  <- seq(0,1,by=0.01)


###################################################
### code chunk number 36: bvpTests.Rnw:493-502
###################################################
xi    <- 0.2
mod1 <- bvpshoot(yini = c(1, NA), yend = c(2, NA), x = x,
                  order = 2, func = prob8, guess = 0)
xi   <- 0.1
mod2 <- bvptwp(yini = c(1, NA), yend = c(2, NA), x = x,
                  order = 2, func = prob8)
xi   <- 0.01
mod3 <- bvptwp(yini = c(1, NA), yend = c(2, NA), x = x,
                  order = 2, func = prob8)


###################################################
### code chunk number 37: prob8
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 8")
curve(2-exp(-1/xi)-exp(-x/xi)/(1-exp(-1/xi)), 
      0, 1, add = TRUE, type = "p")


###################################################
### code chunk number 38: prob8
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 8")
curve(2-exp(-1/xi)-exp(-x/xi)/(1-exp(-1/xi)), 
      0, 1, add = TRUE, type = "p")


###################################################
### code chunk number 39: bvpTests.Rnw:528-531
###################################################
Prob9 <- function(x, y, pars) {
  list(c( y[2], -1/(xi+x^2)*(4*x*y[2]+2*y[1]) ))
}


###################################################
### code chunk number 40: bvpTests.Rnw:533-536
###################################################
xi   <-0.05
mod1 <- bvptwp(yini = c(1/(1+xi), NA), yend = c(1/(1+xi), NA),
               x = seq(-1, 1, by = 0.01), func = Prob9)


###################################################
### code chunk number 41: bvpTests.Rnw:538-541
###################################################
xi   <-0.02
mod2 <- bvptwp(yini = c(1/(1+xi), NA), yend = c(1/(1+xi), NA),
               x = seq(-1, 1, by = 0.01), func = Prob9)


###################################################
### code chunk number 42: bvpTests.Rnw:543-546
###################################################
xi   <-0.01
mod3 <- bvptwp(yini = c(1/(1+xi), NA), yend = c(1/(1+xi), NA),
               x = seq(-1, 1, by = 0.01), func = Prob9)


###################################################
### code chunk number 43: prob9
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 9")
# exact
curve(1/(xi+x^2), -1, 1, type = "p", add = TRUE)


###################################################
### code chunk number 44: prob9
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 9")
# exact
curve(1/(xi+x^2), -1, 1, type = "p", add = TRUE)


###################################################
### code chunk number 45: bvpTests.Rnw:573-576
###################################################
Prob10 <- function(x, y, pars) {
  list( -1/xi*x*y[2] )
}


###################################################
### code chunk number 46: bvpTests.Rnw:578-587
###################################################
xi    <-0.1
mod1 <- bvpshoot(yini = c(0, NA), yend = c(2, NA), 
           order = 2, x = seq(-1, 1, by = 0.01), func=Prob10, guess = 0)
xi   <- 0.05
mod2 <- bvpcol(yini = c(0, NA), yend = c(2, NA), 
           order = 2, x = seq(-1, 1, by = 0.01), func=Prob10)
xi   <- 0.01
mod3 <- bvptwp(yini = c(0, NA), yend = c(2, NA), 
           order = 2, x = seq(-1, 1, by = 0.01), func=Prob10)


###################################################
### code chunk number 47: prob10
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 10")

erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
curve(1+erf(x/sqrt(2*xi))/erf(1/sqrt(2*xi)),
      -1, 1, type = "p", add = TRUE)


###################################################
### code chunk number 48: prob10
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 10")

erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
curve(1+erf(x/sqrt(2*xi))/erf(1/sqrt(2*xi)),
      -1, 1, type = "p", add = TRUE)


###################################################
### code chunk number 49: bvpTests.Rnw:615-618
###################################################
Prob11 <- function(x, y, pars) {
  list(c(y[2], 1/xi * (y[1]-(xi*pi*pi+1)*cos(pi*x)) ))
}


###################################################
### code chunk number 50: bvpTests.Rnw:620-629
###################################################
xi <-0.1
print(system.time(
mod1 <- bvpshoot(yini = c(-1, NA), yend = c(-1, NA), guess = 0,
                   x = seq(-1, 1, by=0.01), func = Prob11, atol = 1e-10)
))
print(system.time(
mod2 <- bvptwp(yini = c(-1, NA), yend = c(-1, NA),
               x = seq(-1, 1, by=0.01), func = Prob11, atol = 1e-10)
))


###################################################
### code chunk number 51: prob11
###################################################
plot(mod1, mod2, which = 1, lty = 1, main = "test problem 11")
curve(cos(pi*x), -1, 1, type = "p", add= TRUE)


###################################################
### code chunk number 52: prob11
###################################################
plot(mod1, mod2, which = 1, lty = 1, main = "test problem 11")
curve(cos(pi*x), -1, 1, type = "p", add= TRUE)


###################################################
### code chunk number 53: bvpTests.Rnw:655-658
###################################################
Prob12 <- function(x, y, pars) {
  list(1/xi * (y[1]-(xi*pi*pi+1)*cos(pi*x)))
}


###################################################
### code chunk number 54: bvpTests.Rnw:660-663
###################################################
xi    <- 0.01
mod1 <- bvpshoot(yini = c(-1, NA), yend = c(0, NA),
           order = 2, x = seq(-1, 1, by = 0.01), func = Prob12, guess = 0)


###################################################
### code chunk number 55: bvpTests.Rnw:665-668
###################################################
xi   <- 0.0025
mod2 <- bvptwp(yini = c(-1, NA), yend = c(0, NA),
           order = 2, x = seq(-1, 1, by = 0.01), func = Prob12)


###################################################
### code chunk number 56: bvpTests.Rnw:670-673
###################################################
xi   <- 0.0001
mod3 <- bvptwp(yini = c(-1, NA), yend = c(0, NA),
           order = 2, x = seq(-1, 1, by = 0.01), func = Prob12)


###################################################
### code chunk number 57: prob12
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 12")
curve(cos(pi*x)+exp((x-1)/sqrt(xi)), 
      -1, 1, type = "p", add = TRUE)


###################################################
### code chunk number 58: prob12
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 12")
curve(cos(pi*x)+exp((x-1)/sqrt(xi)), 
      -1, 1, type = "p", add = TRUE)


###################################################
### code chunk number 59: bvpTests.Rnw:701-704
###################################################
Prob13 <- function(x, y, pars)  {
  list(c( y[2], 1/xi*(y[1]-(xi*pi*pi+1)*cos(pi*x)) ))
}


###################################################
### code chunk number 60: bvpTests.Rnw:706-709
###################################################
xi    <- 0.01
mod1 <- bvpshoot(yini = c(0, NA), yend = c(-1, NA),
                  x = seq(-1, 1, by=0.01), func = Prob13, guess = 0)


###################################################
### code chunk number 61: bvpTests.Rnw:711-714
###################################################
xi   <- 0.0025
mod2 <- bvptwp(yini = c(0, NA), yend = c(-1, NA),
               x = seq(-1, 1, by=0.01), func = Prob13)


###################################################
### code chunk number 62: bvpTests.Rnw:716-719
###################################################
xi    <- 0.0001
mod3 <- bvptwp(yini = c(0, NA), yend = c(-1, NA),
                x = seq(-1, 1, by=0.01), func = Prob13)


###################################################
### code chunk number 63: prob13
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 13")
curve(cos(pi*x)+exp(-(x+1)/sqrt(xi)), 
      -1, 1, type = "p", add = TRUE)


###################################################
### code chunk number 64: prob13
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 13")
curve(cos(pi*x)+exp(-(x+1)/sqrt(xi)), 
      -1, 1, type = "p", add = TRUE)


###################################################
### code chunk number 65: bvpTests.Rnw:746-749
###################################################
Prob14 <- function(x, y, pars)  {
  list(1/xi*(y[1]-(xi*pi*pi+1)*cos(pi*x)))
}


###################################################
### code chunk number 66: bvpTests.Rnw:751-754
###################################################
xi    <- 0.01
mod1 <- bvpshoot(yini = c(0, NA), yend = c(0, NA), 
         order = 2, x = seq(-1, 1, by = 0.01), func = Prob14, guess = 0)


###################################################
### code chunk number 67: bvpTests.Rnw:756-759
###################################################
xi   <- 0.0025
mod2 <- bvptwp(yini = c(0, NA), yend = c(0, NA), 
            order = 2, x = seq(-1, 1, by = 0.01), func = Prob14)


###################################################
### code chunk number 68: bvpTests.Rnw:761-764
###################################################
xi    <- 0.0001
mod3 <- bvptwp(yini = c(0, NA), yend = c(0, NA), 
            order = 2, x = seq(-1, 1, by = 0.01), func = Prob14)


###################################################
### code chunk number 69: prob14
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 14")
curve(cos(pi*x)+exp((x-1)/sqrt(xi))+exp(-(x+1)/sqrt(xi)),
      -1, 1, type = "p", add = TRUE)


###################################################
### code chunk number 70: prob14
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 14")
curve(cos(pi*x)+exp((x-1)/sqrt(xi))+exp(-(x+1)/sqrt(xi)),
      -1, 1, type = "p", add = TRUE)


###################################################
### code chunk number 71: bvpTests.Rnw:787-790
###################################################
Prob15 <- function(x, y, pars)  {
  list(c( y[2], 1/xi*x*y[1] ))
}


###################################################
### code chunk number 72: bvpTests.Rnw:792-797
###################################################
xi <-0.003
print(system.time(
mod1 <- bvpshoot(yini = c(1, NA), yend = c(1, NA),
                   x = seq(-1, 1, by = 0.01), func = Prob15, guess = 0)
))


###################################################
### code chunk number 73: bvpTests.Rnw:799-809
###################################################
xi <- 0.005
print(system.time(
mod2 <- bvptwp(yini = c(1, NA), yend = c(1, NA),
               x = seq(-1, 1, by = 0.01), func = Prob15)
))
xi <- 0.005
print(system.time(
mod3 <- bvpcol(yini = c(1, NA), yend = c(1, NA),
               x = seq(-1, 1, by = 0.01), func = Prob15)
))


###################################################
### code chunk number 74: bvpTests.Rnw:811-816
###################################################
xi <- 0.01
print(system.time(
mod4 <- bvptwp(yini = c(1, NA), yend = c(1, NA),
            x = seq(-1, 1, by = 0.01), func = Prob15)
))


###################################################
### code chunk number 75: prob15
###################################################
plot(mod1, mod2, mod3, mod4, which = 1, lty = 1, main = "test problem 15")


###################################################
### code chunk number 76: prob15
###################################################
plot(mod1, mod2, mod3, mod4, which = 1, lty = 1, main = "test problem 15")


###################################################
### code chunk number 77: bvpTests.Rnw:840-843
###################################################
Prob16 <- function(x, y, pars) {
  list(-1/xi^2*pi^2*y[1]/4 )
}


###################################################
### code chunk number 78: bvpTests.Rnw:845-851
###################################################
xi <-0.11
print(system.time(
mod1 <- bvpshoot(yini = c(0,NA),yend = c(sin(pi/2/xi), NA),
                   x = seq(0, 1, by=0.01), func = Prob16, guess = 0, 
                   order = 2, atol = 1e-10)
))


###################################################
### code chunk number 79: prob16
###################################################
plot(mod1, which = 1, main = "test problem 16", col = "blue")
curve(sin(pi*x/2/xi), 0, 1, type = "p", add = TRUE)


###################################################
### code chunk number 80: prob16
###################################################
plot(mod1, which = 1, main = "test problem 16", col = "blue")
curve(sin(pi*x/2/xi), 0, 1, type = "p", add = TRUE)


###################################################
### code chunk number 81: bvpTests.Rnw:874-879
###################################################
Prob17 <- function(x, y, pars)  {
  list(c( y[2], -3*xi*y[1]/(xi+x^2)^2 ))
}

xseq <- seq(-0.1, 0.1, by = 0.001)


###################################################
### code chunk number 82: bvpTests.Rnw:881-893
###################################################
xi   <- 0.01
mod1 <- bvptwp(yini = c(-0.1/sqrt(xi+0.01), NA),
                yend = c(0.1/sqrt(xi+0.01), NA), x = xseq,
                func = Prob17, atol = 1e-10)
xi   <- 0.001
mod2 <- bvptwp(yini = c(-0.1/sqrt(xi+0.01), NA),
                yend = c(0.1/sqrt(xi+0.01), NA), x = xseq,
                func = Prob17, atol = 1e-8)
xi   <- 0.0001
mod3 <- bvptwp(yini = c(-0.1/sqrt(xi+0.01), NA),
                yend = c(0.1/sqrt(xi+0.01), NA), x = xseq,
                func = Prob17, atol = 1e-8)


###################################################
### code chunk number 83: prob17
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 17")
curve(x/sqrt(xi+x^2), -0.1, 0.1, type = "p", add = TRUE)


###################################################
### code chunk number 84: prob17
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 17")
curve(x/sqrt(xi+x^2), -0.1, 0.1, type = "p", add = TRUE)


###################################################
### code chunk number 85: bvpTests.Rnw:917-922
###################################################
Prob18 <- function(x, y, pars) {
  list( -1/xi*y[2])  
}

xseq<-seq(0,1,by=0.01)


###################################################
### code chunk number 86: bvpTests.Rnw:924-927
###################################################
xi <-0.2
mod1 <- bvpshoot(yini = c(1, NA), yend = c(exp(-1/xi), NA), x = xseq,
              order = 2, func = Prob18, guess = 0, atol = 1e-10)


###################################################
### code chunk number 87: bvpTests.Rnw:929-932
###################################################
xi <- 0.1
mod2 <- bvptwp(yini = c(1, NA), yend = c(exp(-1/xi), NA), x = xseq,
              order = 2, func = Prob18, atol = 1e-10)


###################################################
### code chunk number 88: bvpTests.Rnw:934-937
###################################################
xi <- 0.01
mod3 <- bvptwp(yini = c(1, NA), yend = c(exp(-1/xi), NA), x = xseq,
              order = 2, func = Prob18, atol = 1e-10)


###################################################
### code chunk number 89: prob18
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 18")
curve(exp(-x/xi), 0, 1, type = "p", add = TRUE)


###################################################
### code chunk number 90: prob18
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 18")
curve(exp(-x/xi), 0, 1, type = "p", add = TRUE)


###################################################
### code chunk number 91: bvpTests.Rnw:962-966
###################################################
Prob19 <- function(t, y, pars, ksi) {
  pit = pi*t
  list(c(y[2],(pi/2*sin(pit/2)*exp(2*y[1])-exp(y[1])*y[2])/ksi))
}


###################################################
### code chunk number 92: bvpTests.Rnw:968-979
###################################################
xi    <- 0.05
mod1 <- bvpshoot(yini = c(0, NA), yend = c(0, NA), 
         x = seq(0, 1, by = 0.01), func = Prob19, guess = 0, ksi = xi)
xi <- 0.03
mod2 <- bvptwp(yini = c(0, NA), yend = c(0, NA), 
               x = seq(0, 1, by = 0.01), func = Prob19, ksi = xi, 
               atol = 1e-15)
xi <- 0.005
mod3 <- bvptwp(yini = c(0, NA), yend = c(0, NA), 
                x = seq(0, 1, by = 0.01), func = Prob19, ksi = xi, 
                atol = 1e-10)


###################################################
### code chunk number 93: prob19
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 19")


###################################################
### code chunk number 94: prob19
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 19")


###################################################
### code chunk number 95: bvpTests.Rnw:1002-1005
###################################################
Prob20 <- function(x, y, pars) {
  list( 1/xi *(1-y[2]^2) )
}


###################################################
### code chunk number 96: bvpTests.Rnw:1007-1022
###################################################
xi  <- 0.5
ini <- c(1+xi * log(cosh(0.745/xi)), NA)
end <- c(1+xi * log(cosh(0.255/xi)), NA)
mod1 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by=0.01),
            order = 2, func = Prob20)
xi  <- 0.3
ini <- c(1+xi * log(cosh(0.745/xi)), NA)
end <- c(1+xi * log(cosh(0.255/xi)), NA)
mod2 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by=0.01),
            order = 2, func = Prob20)
xi  <- 0.01
ini <- c(1+xi * log(cosh(0.745/xi)), NA)
end <- c(1+xi * log(cosh(0.255/xi)), NA)
mod3 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by=0.01),
                order = 2, func = Prob20)


###################################################
### code chunk number 97: prob20
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 20")
curve(1+xi * log(cosh((x-0.745)/xi)), 0, 1, add = TRUE, type = "p")


###################################################
### code chunk number 98: prob20
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 20")
curve(1+xi * log(cosh((x-0.745)/xi)), 0, 1, add = TRUE, type = "p")


###################################################
### code chunk number 99: bvpTests.Rnw:1045-1049
###################################################
Prob21 <- function(x, y, pars, xi) {
  list(c( y[2], 1/xi *(y[1]+y[1]^2-exp(-2*x/sqrt(xi))) ))
}
ini <- c(1, NA)


###################################################
### code chunk number 100: bvpTests.Rnw:1051-1063
###################################################
xi  <- 0.2
end <- c(exp(-1/sqrt(xi)), NA)
mod1 <- bvpshoot(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
                   func = Prob21, guess = 0, xi = xi)
xi  <- 0.1
end <- c(exp(-1/sqrt(xi)), NA)
mod2 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
               func = Prob21, xi = xi)
xi  <- 0.01
end <- c(exp(-1/sqrt(xi)), NA)
mod3 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
                func = Prob21, xi = xi)


###################################################
### code chunk number 101: prob21
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 21")
curve(exp(-x/sqrt(xi)), 0, 1, add = TRUE, type = "p")


###################################################
### code chunk number 102: prob21
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 21")
curve(exp(-x/sqrt(xi)), 0, 1, add = TRUE, type = "p")


###################################################
### code chunk number 103: bvpTests.Rnw:1087-1093
###################################################
Prob22 <- function(t, y, pars, xi) {
  list( -1/xi *(y[2]+y[1]^2) )
}

ini <- c(0, NA)
end <- c(1/2, NA)


###################################################
### code chunk number 104: bvpTests.Rnw:1095-1098
###################################################
xi <-0.1
mod1 <- bvpshoot(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
              order = 2, func = Prob22, guess = 0, xi = xi)


###################################################
### code chunk number 105: bvpTests.Rnw:1100-1103
###################################################
xi <-0.05
mod2 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
              order = 2, func = Prob22, xi = xi)


###################################################
### code chunk number 106: bvpTests.Rnw:1105-1108
###################################################
xi <- 0.01
mod3 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
              order = 2, func = Prob22, xi = xi)


###################################################
### code chunk number 107: prob22
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 22")


###################################################
### code chunk number 108: prob22
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 22")


###################################################
### code chunk number 109: bvpTests.Rnw:1131-1137
###################################################
Prob23 <- function(t, y, pars, xi) {
  list(c( y[2], sinh(y[1]/xi)/xi) )
}

ini <- c(0, NA)
end <- c(1, NA)


###################################################
### code chunk number 110: bvpTests.Rnw:1139-1142
###################################################
xi <- 1/5
mod1 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
                func = Prob23, xi = xi)


###################################################
### code chunk number 111: bvpTests.Rnw:1144-1147
###################################################
xi <- 1/7
mod2 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
               func = Prob23, xi = xi)


###################################################
### code chunk number 112: bvpTests.Rnw:1149-1152
###################################################
xi <- 1/9
mod3 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
                func = Prob23, xi = xi)


###################################################
### code chunk number 113: prob23
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 23")


###################################################
### code chunk number 114: prob23
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 23")


###################################################
### code chunk number 115: bvpTests.Rnw:1179-1189
###################################################
Prob24 <- function(t, y, pars, xi) {
  A <- 1+t*t
  AA <- 2*t
  ga <- 1.4
  list((((1+ga)/2 -xi*AA)*y[1]*y[2]-y[2]/y[1]-
       (AA/A)*(1-(ga-1)*y[1]^2/2))/(xi*A*y[1])  )
}

ini <- c(0.9129, NA)
end <- c(0.375, NA)


###################################################
### code chunk number 116: bvpTests.Rnw:1191-1194
###################################################
xi   <- 0.05
mod1 <- bvpshoot(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
            order = 2, func = Prob24, guess = 0.9, xi = xi)


###################################################
### code chunk number 117: bvpTests.Rnw:1196-1200
###################################################
xi <- 0.02
mod2 <- bvpshoot(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
            order = 2,  func = Prob24, guess = 0.9, xi = xi)
attributes(mod2)$roots         # has FAILED: f.root too large!


###################################################
### code chunk number 118: bvpTests.Rnw:1206-1214
###################################################
xi <- 0.02
mod2 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
                order = 2, func = Prob24, xi = xi, 
                xguess = mod2[,1], yguess = t(mod2[,2:3]))
xi <- 0.01
mod3 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
                order = 2, func = Prob24, xi = xi,
                xguess = mod2[,1], yguess = t(mod2[,2:3]))


###################################################
### code chunk number 119: prob24
###################################################
plot(mod1, mod2, mod3,  which = 1, lty = 1, main = "test problem 24")


###################################################
### code chunk number 120: prob24
###################################################
plot(mod1, mod2, mod3,  which = 1, lty = 1, main = "test problem 24")


###################################################
### code chunk number 121: bvpTests.Rnw:1245-1251
###################################################
Prob25 <- function(t, y, pars, xi) {
  list(c( y[2], -1/xi *(y[1]*y[2]-y[1]) ))
}

ini <- c(-1/3 ,NA)
end <- c(1/3, NA)


###################################################
### code chunk number 122: bvpTests.Rnw:1253-1262
###################################################
xi   <- 0.1
mod1 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
               func = Prob25, xi = xi)
xi   <- 0.01
mod2 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
               func = Prob25, xi = xi)
xi   <- 0.001
mod3 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
               func = Prob25, xi = xi)


###################################################
### code chunk number 123: prob25
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 25")


###################################################
### code chunk number 124: prob25
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 25")


###################################################
### code chunk number 125: bvpTests.Rnw:1286-1292
###################################################
Prob26 <- function(t, y, pars, xi)  {
  list(  -1/xi *(y[1]*y[2]-y[1]) )
}

ini <- c(1, NA)
end <- c(-1/3, NA)


###################################################
### code chunk number 126: bvpTests.Rnw:1294-1297
###################################################
xi   <- 0.1
mod1 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
               order = 2, func = Prob26, xi = xi)


###################################################
### code chunk number 127: bvpTests.Rnw:1299-1302
###################################################
xi   <- 0.02
mod2 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
               order = 2, func = Prob26, xi = xi)


###################################################
### code chunk number 128: bvpTests.Rnw:1304-1307
###################################################
xi   <- 0.005
mod3 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
               order = 2, func = Prob26, xi = xi)


###################################################
### code chunk number 129: prob26
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 26")


###################################################
### code chunk number 130: prob26
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 26")


###################################################
### code chunk number 131: bvpTests.Rnw:1329-1335
###################################################
Prob27 <- function(t, y, pars, xi) {
  list(c( y[2], -1/xi *(y[1]*y[2]-y[1]) ))
}

ini <- c(1, NA)
end <- c(1/3, NA)


###################################################
### code chunk number 132: bvpTests.Rnw:1337-1340
###################################################
xi   <- 0.1
mod1 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
               func = Prob27, xi = xi)


###################################################
### code chunk number 133: bvpTests.Rnw:1342-1345
###################################################
xi   <- 0.02
mod2 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
               func = Prob27, xi = xi)


###################################################
### code chunk number 134: bvpTests.Rnw:1347-1350
###################################################
xi   <- 0.005
mod3 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
               func = Prob27, xi = xi)


###################################################
### code chunk number 135: prob27
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 27")


###################################################
### code chunk number 136: prob27
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 27")


###################################################
### code chunk number 137: bvpTests.Rnw:1373-1379
###################################################
Prob28 <- function(t, y, pars, xi)  {
  list(  -1/xi *(y[1]*y[2]-y[1]))
}

ini <- c(1, NA)
end <- c(3/2, NA)


###################################################
### code chunk number 138: bvpTests.Rnw:1381-1384
###################################################
xi   <- 0.1
mod1 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
               order = 2, func = Prob28, xi = xi)


###################################################
### code chunk number 139: bvpTests.Rnw:1386-1389
###################################################
xi   <-0.02
mod2 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
               order = 2, func = Prob28, xi = xi)


###################################################
### code chunk number 140: bvpTests.Rnw:1391-1394
###################################################
xi <-0.005
mod3 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
               order = 2, func = Prob28, xi = xi)


###################################################
### code chunk number 141: prob28
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 28")


###################################################
### code chunk number 142: prob28
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 28")


###################################################
### code chunk number 143: bvpTests.Rnw:1416-1422
###################################################
Prob29 <- function(t, y, pars, xi)  {
  list(c( y[2], -1/xi *(y[1]*y[2]-y[1]) ))
}

ini <- c(0,NA)
end <- c(3/2,NA)


###################################################
### code chunk number 144: bvpTests.Rnw:1424-1427
###################################################
xi   <- 0.1
mod1 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
               func = Prob29, xi = xi)


###################################################
### code chunk number 145: bvpTests.Rnw:1429-1432
###################################################
xi   <- 0.02
mod2 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
               func = Prob29, xi = xi)


###################################################
### code chunk number 146: bvpTests.Rnw:1434-1437
###################################################
xi   <- 0.005
mod3 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
               func = Prob29, xi = xi)


###################################################
### code chunk number 147: prob29
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 29")


###################################################
### code chunk number 148: prob29
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 29")


###################################################
### code chunk number 149: bvpTests.Rnw:1461-1466
###################################################
Prob30 <- function(t, y, pars, xi)  {
  list( -1/xi *(y[1]*y[2]-y[1]) )
}
ini <- c(-7/6, NA)
end <- c(3/2, NA)


###################################################
### code chunk number 150: bvpTests.Rnw:1468-1477
###################################################
xi   <- 0.1
mod1 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
               order = 2, func = Prob30, xi = xi)
xi   <- 0.02
mod2 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
               order = 2, func = Prob30, xi = xi)
xi   <- 0.01
mod3 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
               order = 2, func = Prob30, xi = xi)


###################################################
### code chunk number 151: prob30
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 30")


###################################################
### code chunk number 152: prob30
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 30")


###################################################
### code chunk number 153: bvpTests.Rnw:1505-1518
###################################################
Prob31 <- function(t, Y, pars)  {
  with (as.list(Y), {
    dy    <- sin(Tet)
    dTet  <- M
    dM    <- -Q/xi
    T     <- 1/cos (Tet) +xi*Q*tan(Tet)
    dQ    <- 1/xi*((y-1)*cos(Tet)-M*T)
    list(c( dy, dTet, dM, dQ))
  })
}

ini <- c(y = 0, Tet = NA, M = 0, Q = NA)
end <- c(y = 0, Tet = NA, M = 0, Q = NA)


###################################################
### code chunk number 154: bvpTests.Rnw:1522-1531
###################################################
xi <-0.1
mod1 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
               func = Prob31, atol = 1e-10)
xi <- 0.05
mod2 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
               func = Prob31, atol = 1e-10)
xi <- 0.01
mod3 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
               func = Prob31, atol = 1e-10)


###################################################
### code chunk number 155: prob31
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 31")


###################################################
### code chunk number 156: prob31
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 31")


###################################################
### code chunk number 157: bvpTests.Rnw:1555-1560
###################################################
Prob32 <- function(t, y, pars, xi) {
  list(1/xi*(y[2]*y[3]-y[1]*y[4]))
}
ini <- c(0, 0, NA, NA)
end <- c(1, 0, NA, NA)


###################################################
### code chunk number 158: bvpTests.Rnw:1562-1571
###################################################
xi  <- 0.01
mod1 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
              order = 4, func = Prob32, xi = xi)
xi   <- 0.002
mod2 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
              order = 4, func = Prob32, xi = xi)
xi   <- 0.0001
mod3 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
              order = 4, func = Prob32, xi = xi)


###################################################
### code chunk number 159: prob32
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 32")


###################################################
### code chunk number 160: prob32
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 32")


###################################################
### code chunk number 161: bvpTests.Rnw:1597-1603
###################################################
Prob33 <- function(t, z, pars, xi) {
  list(c( z[2], z[3], z[4], 1/xi*(z[1]*z[4]-z[5]*z[6]),
          z[6], 1/xi*(z[5]*z[2]-z[1]*z[6])))
}
ini <- c(0, 0, NA, NA, -1, NA)
end <- c(0, 0, NA, NA,  1, NA)


###################################################
### code chunk number 162: bvpTests.Rnw:1605-1614
###################################################
xi  <- 0.1
mod1 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
              func = Prob33, xi = xi)
xi  <- 0.01
mod2 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
              func = Prob33, xi = xi)
xi  <- 0.001
mod3 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
              func = Prob33, xi = xi)


###################################################
### code chunk number 163: prob33
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 33")


###################################################
### code chunk number 164: prob33
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 33")


###################################################
### code chunk number 165: bvpTests.Rnw:1640-1645
###################################################
Prob34 <- function(t, y, pars, xi) {
  list(-xi*exp(y[1]))
}
ini <- c(0, NA)
end <- c(0, NA)


###################################################
### code chunk number 166: bvpTests.Rnw:1647-1650
###################################################
xi  <- 0.1
mod1 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
              order = 2, func = Prob34, xi = xi)


###################################################
### code chunk number 167: bvpTests.Rnw:1652-1655
###################################################
xi   <- 0.01
mod2 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
              order = 2, func = Prob34, xi = xi)


###################################################
### code chunk number 168: bvpTests.Rnw:1657-1660
###################################################
xi   <- 0.001
mod3 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
              order = 2, func = Prob34, xi = xi)


###################################################
### code chunk number 169: prob34
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 34")


###################################################
### code chunk number 170: prob34
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 34")


###################################################
### code chunk number 171: bvpTests.Rnw:1685-1690
###################################################
Prob35 <- function(x, y, pars, xi) {
  list(c( y[2], 1/xi*(x * y[2]-y[1])))
}
ini <- c(1, NA)
end <- c(2, NA)


###################################################
### code chunk number 172: bvpTests.Rnw:1692-1695
###################################################
xi  <- 1
mod1 <- bvptwp(yini = ini, yend = end, x = seq(-1, 1, by = 0.05),
              func = Prob35, xi = xi)


###################################################
### code chunk number 173: bvpTests.Rnw:1697-1700
###################################################
xi   <- 0.1
mod2 <- bvptwp(yini = ini, yend = end, x = seq(-1, 1, by = 0.05),
              func = Prob35, xi = xi)


###################################################
### code chunk number 174: bvpTests.Rnw:1702-1705
###################################################
xi   <- 0.01
mod3 <- bvptwp(yini = ini, yend = end, x = seq(-1, 1, by = 0.05),
              func = Prob35, xi = xi)


###################################################
### code chunk number 175: prob35
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 35")


###################################################
### code chunk number 176: prob35
###################################################
plot(mod1, mod2, mod3, which = 1, lty = 1, main = "test problem 35")


