### R code from vignette source 'bvpSolve.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: preliminaries
###################################################
library("bvpSolve")
options(prompt = " ")
options(continue = "  ")
options(width=70)


###################################################
### code chunk number 2: bvpSolve.Rnw:283-288
###################################################
fun<- function(x,y,pars) {
 list(c(y[2],
   1/ks * (-x*y[2]+y[1]-(1+ks*pi*pi)*cos(pi*x)-pi*x*sin(pi*x)))
     )
}


###################################################
### code chunk number 3: bvpSolve.Rnw:292-306
###################################################
ks <- 0.1
x  <- seq(-1, 1, by = 0.01)
print(system.time(
sol1  <- bvpshoot(yini = c(-1, NA), yend = c(1, NA),
                  x = x, func = fun, guess = 0)
))

print(system.time(
sol2  <- bvptwp(yini = c(-1, NA), yend = c(1, NA), x = x, func = fun)
))

print(system.time(
sol3  <- bvpcol(yini = c(-1, NA), yend = c(1, NA), x = x, func = fun)
))


###################################################
### code chunk number 4: pr7a
###################################################
plot(sol2[,1], sol2[,3], type = "l", main = "test problem 7, ksi=0.1",
     lwd = 2, col = "red")
points(sol1[,1], sol1[,3], col = "green", pch = "x")
legend("topright", c("bvptwp", "bvpshoot"),
        lty = c(1, NA, NA), pch = c(NA, 1, 3), col = c("red", "green"))


###################################################
### code chunk number 5: pr7a
###################################################
plot(sol2[,1], sol2[,3], type = "l", main = "test problem 7, ksi=0.1",
     lwd = 2, col = "red")
points(sol1[,1], sol1[,3], col = "green", pch = "x")
legend("topright", c("bvptwp", "bvpshoot"),
        lty = c(1, NA, NA), pch = c(NA, 1, 3), col = c("red", "green"))


###################################################
### code chunk number 6: bvpSolve.Rnw:345-355
###################################################
ks <- 0.0005

print(system.time(
sol3  <- bvptwp(yini = c(-1, NA), yend = c(1, NA), x = seq(-1, 1, by = 0.01),
           func = fun, xguess = sol2[,1], yguess = t(sol2[,-1]))
))
print(system.time(
sol3b <- bvpcol(yini = c(-1, NA), yend = c(1, NA), x = seq(-1, 1, by = 0.01),
           func = fun, xguess = sol2[,1], yguess = t(sol2[,-1]))
))


###################################################
### code chunk number 7: bvpSolve.Rnw:359-364
###################################################
ks <- 0.0001
print(system.time(
sol3c <- bvpcol(yini = c(-1, NA), yend = c(1, NA), x = seq(-1, 1, by = 0.01),
           func = fun, yguess = sol3b)
))


###################################################
### code chunk number 8: pr7b
###################################################
plot(sol3, type = "l", lwd = 2, col = "red")



###################################################
### code chunk number 9: pr7b
###################################################
plot(sol3, type = "l", lwd = 2, col = "red")



###################################################
### code chunk number 10: bvpSolve.Rnw:415-425
###################################################
fsub <- function (t,Y,pars)  { 
  return(list(c(f1 = Y[2],
                f2 = (Y[1]*Y[4] - Y[3]*Y[2])/eps,
                f3 = Y[4],
                f4 = Y[5],
                f5 = Y[6],
                f6 = (-Y[3]*Y[6] - Y[1]*Y[2])/eps)))
}
eps <- 0.001
x <- seq(0, 1, len = 100)


###################################################
### code chunk number 11: bvpSolve.Rnw:429-439
###################################################
print(system.time(
Soltwp <- bvptwp(x = x, func = fsub,
           yini = c(y1 = -1, y2 = NA, y3 = 0, y4 = 0, y5 = NA, y6 = NA),
           yend = c(1, NA, 0, 0, NA, NA))
))
print(system.time(
Solcol <- bvpcol(x = x, func = fsub,
           yini = c(y1 = -1, y2 = NA, y3 = 0, y4 = 0, y5 = NA, y6 = NA),
           yend = c(1, NA, 0, 0, NA, NA))
))


###################################################
### code chunk number 12: bvpSolve.Rnw:447-449
###################################################
diagnostics(Soltwp)
diagnostics(Solcol)


###################################################
### code chunk number 13: swirl
###################################################
pairs(Soltwp, main = "swirling flow III, eps=0.01", col = "blue")


###################################################
### code chunk number 14: swirl
###################################################
pairs(Soltwp, main = "swirling flow III, eps=0.01", col = "blue")


###################################################
### code chunk number 15: bvpSolve.Rnw:495-501
###################################################
fsubhigh <- function (t,Y,pars)  { 
  return(list(c(d2g = (Y[1]*Y[4] - Y[3]*Y[2])/eps,
                d4f = (-Y[3]*Y[6] - Y[1]*Y[2])/eps)))
}
eps <- 0.001
x   <- seq(0, 1, len = 100)


###################################################
### code chunk number 16: bvpSolve.Rnw:505-510
###################################################
print(system.time(
Solcol2 <- bvpcol(x = x, func = fsubhigh, order = c(2, 4),
           yini = c(y1 = -1, y2 = NA, y3 = 0, y4 = 0, y5 = NA, y6 = NA),
           yend = c(1, NA, 0, 0, NA, NA))
))


###################################################
### code chunk number 17: bvpSolve.Rnw:515-517
###################################################
max(abs(Solcol2-Solcol))
diagnostics(Solcol2)


###################################################
### code chunk number 18: bvpSolve.Rnw:522-527
###################################################
print(system.time(
Soltwp2 <- bvptwp(x = x, func = fsubhigh, order = c(2, 4),
           yini = c(y1 = -1, y2 = NA, y3 = 0, y4 = 0, y5 = NA, y6 = NA),
           yend = c(1, NA, 0, 0, NA, NA))
))


###################################################
### code chunk number 19: bvpSolve.Rnw:531-532
###################################################
max(abs(Soltwp2- Soltwp))


###################################################
### code chunk number 20: bvpSolve.Rnw:543-546
###################################################
eps <- 0.0001
xguess <- Soltwp[,1]
yguess <- t(Soltwp[,2:7])


###################################################
### code chunk number 21: bvpSolve.Rnw:548-553
###################################################
print(system.time(
  Sol2 <- bvptwp(x = x, func = fsub, xguess = xguess, yguess = yguess, 
          yini = c(y1 = -1, y2 = NA, y3 = 0, y4 = 0, y5 = NA, y6 = NA),
          yend = c(1,            NA,      0,      0,      NA,      NA))
))


###################################################
### code chunk number 22: bvpSolve.Rnw:558-563
###################################################
print(system.time(
  Sol2b <- bvpcol(x = x, func = fsub, yguess = Solcol, 
          yini = c(y1 = -1, y2 = NA, y3 = 0, y4 = 0, y5 = NA, y6 = NA),
          yend = c(1,            NA,      0,      0,      NA,      NA))
))


###################################################
### code chunk number 23: swirl2
###################################################
plot(Sol2, col = "darkred", type = "l", lwd = 2)
mtext(outer = TRUE, side = 3, line = -1.5, cex = 1.5,
   "swirling flow III, eps=0.0001")


###################################################
### code chunk number 24: swirl2
###################################################
plot(Sol2, col = "darkred", type = "l", lwd = 2)
mtext(outer = TRUE, side = 3, line = -1.5, cex = 1.5,
   "swirling flow III, eps=0.0001")


###################################################
### code chunk number 25: swirl3
###################################################
plot (Solcol, Sol2b, which = c("y1","y4"), lty = 1)


###################################################
### code chunk number 26: swirl3
###################################################
plot (Solcol, Sol2b, which = c("y1","y4"), lty = 1)


###################################################
### code chunk number 27: bvpSolve.Rnw:610-614
###################################################
obsdat <- matrix (ncol = 2, data = 
      c(seq(0, 1, 0.2), c(-1, -0.4, -0.1, 0.1, 0.4, 1)))
colnames (obsdat) <- c("x", "y1")
obsdat


###################################################
### code chunk number 28: swirl4
###################################################
plot(Solcol, Sol2b, obs = obsdat)


###################################################
### code chunk number 29: swirl4
###################################################
plot(Solcol, Sol2b, obs = obsdat)


###################################################
### code chunk number 30: bvpSolve.Rnw:655-665
###################################################
musn <- function(x,Y,pars) {
  with (as.list(Y),  {
    du <- 0.5 * u * (w - u) /v
    dv <- -0.5 * (w - u)
    dw <- (0.9 - 1000 * (w - y) - 0.5 * w * (w - u)) /z
    dz <- 0.5 * (w - u)
    dy <- -100 * (y - w)
    return(list(c(du, dv, dw, dz, dy)))
  })
}


###################################################
### code chunk number 31: bvpSolve.Rnw:676-677
###################################################
init <- c(u = 1, v = 1, w = 1, z = -10, y = NA)


###################################################
### code chunk number 32: bvpSolve.Rnw:688-689
###################################################
yend  <- function (Y, yini, pars)  with (as.list(Y), w-y)


###################################################
### code chunk number 33: bvpSolve.Rnw:698-702
###################################################
print(system.time(
sol   <-bvpshoot(yini = init, x = seq(0, 1, by = 0.05), func = musn,
           yend = yend, guess = 1, atol = 1e-10, rtol = 0)
))


###################################################
### code chunk number 34: musn
###################################################
plot(sol, type = "l", lwd = 2)
mtext(outer = TRUE, side = 3, line = -1.5, cex = 1.5, "musn")


###################################################
### code chunk number 35: musn
###################################################
plot(sol, type = "l", lwd = 2)
mtext(outer = TRUE, side = 3, line = -1.5, cex = 1.5, "musn")


###################################################
### code chunk number 36: bvpSolve.Rnw:725-734
###################################################
bound <- function(i,y,pars) {
  with (as.list(y), {
    if (i ==1) return (u-1)
    if (i ==2) return (v-1)
    if (i ==3) return (w-1)
    if (i ==4) return (z+10)
    if (i ==5) return (w-y)
 })
}


###################################################
### code chunk number 37: bvpSolve.Rnw:739-744
###################################################
xguess <- seq(0, 1, len = 5)
yguess <- matrix(ncol = 5, (rep(c(1, 1, 1, -10, 0.91), times = 5)) )
rownames(yguess) <- c("u","v","w","z","y")
xguess
yguess


###################################################
### code chunk number 38: bvpSolve.Rnw:752-764
###################################################
print(system.time(
Sol <- bvptwp(yini = NULL, x = x, func = musn, bound = bound,
              xguess = xguess, yguess = yguess, leftbc = 4,
              atol = 1e-10)
))

print(system.time(
Sol2 <- bvpcol(yini = NULL, x = x, func = musn, bound = bound,
              xguess = xguess, yguess = yguess, leftbc = 4,
              atol = 1e-10)
))



###################################################
### code chunk number 39: bvpSolve.Rnw:789-792
###################################################
mathieu <- function(x,y,lambda)
    list(c(y[2],
        -(lambda - 10 * cos(2 * x)) * y[1]))


###################################################
### code chunk number 40: bvpSolve.Rnw:797-800
###################################################
init <- c(1, 0)
sol  <- bvpshoot(yini = init, yend = c(NA, 0), x = seq(0, pi, by = 0.01),
        func = mathieu, extra = 15)


###################################################
### code chunk number 41: mat
###################################################
plot(sol[,1:2])
mtext(outer = TRUE, side = 3, line = -1.5, cex = 1.5, "mathieu")


###################################################
### code chunk number 42: mat
###################################################
plot(sol[,1:2])
mtext(outer = TRUE, side = 3, line = -1.5, cex = 1.5, "mathieu")


###################################################
### code chunk number 43: bvpSolve.Rnw:819-820
###################################################
attr(sol, "roots")  # root gives the value of "lambda" (17.10683)


###################################################
### code chunk number 44: bvpSolve.Rnw:834-838
###################################################
mathieu2 <- function(x,y,p)
    list(c(y[2],
        -(y[3] - 10 * cos(2 * x)) * y[1],
        0) )


###################################################
### code chunk number 45: bvpSolve.Rnw:848-851
###################################################
Sol <- bvptwp (yini = c(y = 1, dy = 0, lambda = NA), yend = c(NA, 0, NA),
        x = seq(0, pi, by = 0.01), func = mathieu2, xguess = c(0, pi),
        yguess = matrix(nrow = 3, data = rep(15, 6)) )


###################################################
### code chunk number 46: mat2
###################################################
plot(Sol, type = "l", lwd = 2)
mtext(outer = TRUE, side = 3, line = -1.5, cex = 1.5, 
     "mathieu - solved using bvptwp")


###################################################
### code chunk number 47: mat2
###################################################
plot(Sol, type = "l", lwd = 2)
mtext(outer = TRUE, side = 3, line = -1.5, cex = 1.5, 
     "mathieu - solved using bvptwp")


###################################################
### code chunk number 48: bvpSolve.Rnw:892-895
###################################################
nerve <- function (t,y,T)
  list(c( 3 * T * (y[1] + y[2] - 1/3 * (y[1]^3) - 1.3),
        (-1/3) * T * (y[1] - 0.7 + 0.8 * y[2]) ))


###################################################
### code chunk number 49: bvpSolve.Rnw:898-902
###################################################
res<- function (Y,yini,T)
  c(Y[1] - yini[1],
    Y[2] - yini[2],
    T*(-1/3) * (yini[1] - 0.7 + 0.8 * yini[2]) - 1)


###################################################
### code chunk number 50: bvpSolve.Rnw:907-910
###################################################
yini <- c(y1 = NA, y2 = NA)
sol  <- bvpshoot(yini = yini, x = seq(0, 1, by = 0.01),
        func = nerve, guess = c(0.5,0.5), yend = res, extra = 2 * pi)


###################################################
### code chunk number 51: bvpSolve.Rnw:913-914
###################################################
attributes(sol)$root


###################################################
### code chunk number 52: bvpSolve.Rnw:929-936
###################################################
nerve3 <- function (t,y,p)
  list(c( 3 * y[3] * (y[1] + y[2] - 1/3 * (y[1]^3) - 1.3),
        (-1/3) * y[3] * (y[1] - 0.7 + 0.8 * y[2]) ,
        0,
        0,
        0)
  )


###################################################
### code chunk number 53: bvpSolve.Rnw:940-947
###################################################
bound <- function(i,y,p) {
 if (i == 1) return ( y[3]*(-1/3) * (y[1] - 0.7 + 0.8 * y[2]) - 1 )
 if (i == 2) return ( y[1] - y[4] )           
 if (i == 3) return ( y[2] - y[5] )           # left bnd
 if (i == 4) return ( y[1] - y[4] )           # right bnd
 if (i == 5) return ( y[2] - y[5] )
}


###################################################
### code chunk number 54: bvpSolve.Rnw:955-960
###################################################
xguess = seq(0, 1, by = 0.1)
yguess = matrix(nrow = 5, ncol = length(xguess), data = 5.)
yguess[1,] <- sin(2 * pi * xguess)
yguess[2,] <- cos(2 * pi * xguess)
rownames(yguess) <- c("y1", "y2", "T", "y1ini", "y2ini")


###################################################
### code chunk number 55: bvpSolve.Rnw:964-967
###################################################
Sol  <- bvptwp(func = nerve3, bound = bound, x = seq(0, 1, by = 0.01), 
        ynames = c("y", "dy", "T", "yi", "yj"),
        leftbc = 3, xguess = xguess, yguess = yguess)


###################################################
### code chunk number 56: bvpSolve.Rnw:970-971
###################################################
Sol[1,]


###################################################
### code chunk number 57: nerve
###################################################
plot(Sol, type = "l", lwd = 2, col = "darkblue")


###################################################
### code chunk number 58: nerve
###################################################
plot(Sol, type = "l", lwd = 2, col = "darkblue")


###################################################
### code chunk number 59: bvpSolve.Rnw:1009-1032
###################################################
fluid<-function(t, y, pars, R) {
 P    <- 0.7*R
 with(as.list(y),   {
  df = f1                       #f'
  df1= f2                       #f''
  df2= R * (f1^2 - f*f2)-A      #f'''
  dh = h1
  dh1= -R * f * h1 - 1
  dO = O1
  dO1= -P * f * O1             
  dA = 0                    # the constant to be estimated
  return(list(c(df, df1, df2, dh, dh1, dO, dO1, dA)))
 })
}

times  <- seq(0, 1, by = 0.01)
yini   <- c(f = 0, f1 = 0, f2 = NA, h = 0, h1 = NA, O = 0, O1 = NA, A = NA)
yend   <- c(1,          0,      NA,     0,      NA,     1,      NA,     NA)

print (system.time(
Solcol1 <- bvpcol(func=fluid, x=times, parms=NULL, R=10000,
                 yini = yini, yend=yend)
))                 


###################################################
### code chunk number 60: bvpSolve.Rnw:1039-1058
###################################################
fluidHigh <- function(t, y, pars, R) {
 P    <- 0.7 * R
 with(as.list(y),   {
  d3f = R * (f1^2 - f * f2) - A      #f'''
  d2h = -R * f * h1 - 1
  d2O = -P * f * O1               
  dA  = 0
  return(list(c(d3f, d2h, d2O, dA)))
 })
}

times  <- seq(0, 1, by = 0.01)
yini   <- c(f = 0, f1 = 0, f2 = NA, h = 0, h1 = NA, O = 0, O1 = NA, A = NA)
yend   <- c(1,          0,      NA,     0,      NA,     1,      NA,     NA)

print (system.time(
Solcol2 <- bvpcol(func = fluidHigh, x = times, parms = NULL, R = 10000, 
                  order = c(3, 2, 2, 1), yini = yini, yend=yend)
))


###################################################
### code chunk number 61: bvpSolve.Rnw:1064-1066
###################################################
head(Solcol1, n = 3)
head(Solcol2, n = 3)


###################################################
### code chunk number 62: fluid
###################################################
plot(Solcol1, main="Fluid injection problem", 
     which = "f1", type = "l", lwd = 2)


###################################################
### code chunk number 63: fluid
###################################################
plot(Solcol1, main="Fluid injection problem", 
     which = "f1", type = "l", lwd = 2)


###################################################
### code chunk number 64: bvpSolve.Rnw:1112-1125
###################################################
multip <- function (x, y, p) {
  list(c((y[2] - 1)/2, 
         (y[1]*y[2] - x)/mu))
}

bound <- function (i, y, p) {
  if (i == 1) y[2] -1         # at x=0.5: y2=1
  else y[1]                   # at x=  1: y1=0
}

mu  <- 0.1
sol <- bvpcol(func = multip, bound = bound, 
              x = seq(0, 1, 0.01), posbound = c(0.5, 1))


###################################################
### code chunk number 65: bvpSolve.Rnw:1128-1129
###################################################
sol[sol[,1] %in% c(0.5,1),]


###################################################
### code chunk number 66: multipoint
###################################################
plot(sol)


###################################################
### code chunk number 67: multipoint
###################################################
plot(sol)


###################################################
### code chunk number 68: bvpSolve.Rnw:1172-1185
###################################################
Elastica <- function (x, y, pars) {

  list( c(cos(y[3]),
          sin(y[3]),
          y[4],
          y[5] * cos(y[3]),
          0))
}

Sol <- bvpcol(func = Elastica,
              yini = c(x = 0,  y = 0, p = NA,   k = 0,  F = NA),
              yend = c(x = NA, y = 0, p = -pi/2,k = NA, F = NA),
              x = seq(0, 0.5, len = 16))


###################################################
### code chunk number 69: elastica
###################################################
plot(Sol)


###################################################
### code chunk number 70: elastica
###################################################
plot(Sol)


###################################################
### code chunk number 71: bvpSolve.Rnw:1213-1223
###################################################
jacfunc <- function (x, y, pars) {
      Jac <- matrix(nrow = 5, ncol = 5, data = 0)
      Jac[3,4] <- 1.0
      Jac[4,4] <- 1.0
      Jac[1,3] <- -sin(y[3])
      Jac[2,3] <- cos(y[3])
      Jac[4,3] <- -y[5] * sin(y[3])
      Jac[4,5] <- Jac[2,3]
      Jac
}


###################################################
### code chunk number 72: bvpSolve.Rnw:1226-1232
###################################################
bound <- function (i, y, pars)  {
    if      (i <= 2) return(y[i])
    else if (i == 3) return(y[4])
    else if (i == 4) return(y[2])
    else if (i == 5) return(y[3] + pi/2)
}


###################################################
### code chunk number 73: bvpSolve.Rnw:1234-1242
###################################################
jacbound <- function(i, y, pars)  {
    JJ <- rep(0, 5)
         if (i <= 2) JJ[i] =1.0
    else if (i == 3) JJ[4] =1.0
    else if (i == 4) JJ[2] =1.0
    else if (i == 5) JJ[3] =1.0
    JJ
}


###################################################
### code chunk number 74: bvpSolve.Rnw:1247-1251
###################################################
Sol4 <- bvpcol(leftbc = 3, ynames = c("x", "y", "p", "k", "F"),
              func = Elastica, jacfunc = jacfunc,
              bound = bound, jacbound = jacbound,
              x = seq(0, 0.5, len=16))


###################################################
### code chunk number 75: bvpSolve.Rnw:1480-1484
###################################################
outF <- bvpcol(ncomp = 5,
               x = seq(0, 0.5, len = 16), leftbc = 3, func = "fsub", 
               jacfunc = "dfsub", bound = "gsub", jacbound = "dgsub",
               dllname = "bvpSolve")


###################################################
### code chunk number 76: bvpSolve.Rnw:1543-1547
###################################################
fun <- function(t,y,pars)
  list(c( y[2],
        - a * p * y[1]/(p + t*t)^2
        ))


###################################################
### code chunk number 77: bvpSolve.Rnw:1550-1552
###################################################
p    <- 1e-5
a    <- 3


###################################################
### code chunk number 78: bvpSolve.Rnw:1557-1561
###################################################
sol  <- bvptwp(yini = c(y = -0.1/sqrt(p+0.01), dy = NA),
               yend = c(     0.1/sqrt(p+0.01),      NA),
               x = seq(-0.1, 0.1, by = 0.001),
               func = fun)


###################################################
### code chunk number 79: linear
###################################################
plot(sol, type = "l")


###################################################
### code chunk number 80: linear
###################################################
plot(sol, type = "l")


###################################################
### code chunk number 81: bvpSolve.Rnw:1649-1650
###################################################
parms <- c(a = 3, p = 1e-7)


###################################################
### code chunk number 82: bvpSolve.Rnw:1669-1681
###################################################
Out  <- NULL
x    <- seq(-0.1, 0.1, by = 0.001)
pseq <- 10^-seq(0, 6, 0.5)

for (pp in pseq) {
  parms[2] <- pp
  outFor <- bvptwp(ncomp = 2, x = x, leftbc = 1, 
     initfunc = "initbnd", parms = parms, func = "funbnd", 
     jacfunc = "dfbnd", bound = "gbnd", jacbound = "dgbnd",
     allpoints = FALSE, dllname = "bvpSolve")
  Out <- cbind(Out, outFor[,2])
}


###################################################
### code chunk number 83: bvpSolve.Rnw:1686-1686
###################################################



###################################################
### code chunk number 84: linf
###################################################
matplot(x, Out, type = "l")
legend("topleft", legend = log10(pseq), title = "logp",
  col = 1 : length(pseq), lty = 1 : length(pseq), cex = 0.6)


###################################################
### code chunk number 85: linf
###################################################
matplot(x, Out, type = "l")
legend("topleft", legend = log10(pseq), title = "logp",
  col = 1 : length(pseq), lty = 1 : length(pseq), cex = 0.6)


