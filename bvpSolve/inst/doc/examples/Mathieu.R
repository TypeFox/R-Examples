## =============================================================================
## Find the 4th eigenvalue of Mathieu's equation:
## y''+(lam-10cos2t)y=0   on the interval [0,pi]
## y(0)=1, y'(0)=0  and y'(pi)=0
##
## 2nd order problem is rewritten as:
## dy=y2
## dy2= -(lam-10cos(2t))*y
## =============================================================================

require(bvpSolve)
require(rootSolve)

mathieu<- function(t, y, lambda = 15)  {
 list(c(y[2], -(lambda-10*cos(2*t)) * y[1]))
}

## =============================================================================
## Solution method 1:  shooting  
## =============================================================================
x = seq(0, pi, by = 0.01)

init <- c(1, 0)
sol  <- bvpshoot(yini = init, yend = c(NA, 0), x = x,
          func = mathieu, guess = NULL,  extra = 15)
plot(sol[,1:2])

## =============================================================================
## Solution method 2: multiroot + bvptwp
## =============================================================================

cost <- function(X)
{  sol<- bvptwp(yini = c(1, NA), yend = c(NA, 0), 
          x = c(0, pi), parms = X,
          func = mathieu)
  return(sol[2,3])  # y2[0]=0
}

# find the root
lam <- multiroot(f = cost, start = 15)

# solve the mode with this root...
Sol<- bvptwp(yini = c(1,NA), yend = c(NA, 0),
        x = x, parms = lam$root,
        func = mathieu, atol = 1e-10)
lines(Sol, col = "red")

## =============================================================================
## Solution method 3 augmented equations...
## =============================================================================

mathieu2<- function(t,y,parms)
{
 list(c(y[2],
        -(y[3]-10*cos(2*t))*y[1],
        0 ))
}

### initial guess = 1 ####

# Solvable with bvpshoot
init <- c(y = 1,dy = 0, lambda = NA)
sol1 <- bvpshoot(yini = init, yend = c(NA, 0, NA), x = x,
          func = mathieu2, guess = 1)
plot(sol1)

jac <- function(x, y ,p) {
  df <- matrix(nr = 3, nc = 3, 0)
  df[1,2] <- 1
  df[2,1] <- -(y[3]-10*cos(2*x))
  df[2,3] <- -y[1]
  df
}
xguess <-  c(0, 1, 2*pi)
yguess <- matrix(nr = 3, rep(1, 9))
rownames(yguess) <- c("y", "dy", "lambda")

# only works for bvptwp if yguess is not 0!...
print(system.time(
sol1b <- bvptwp(yini = init, yend = c(NA, 0, NA), x = x,
        func = mathieu2, jacfunc = jac, xguess = xguess,
        yguess = yguess)
))

plot(sol1b, type="l",lwd=2)

### initial guess = 17 ####
xguess <-  c(0,1,2*pi)
yguess <- matrix(nr=3,rep(17,9))

print(system.time(
sol2 <- bvpshoot(yini = init, yend = c(NA, 0, NA), x = x,
        func = mathieu2, jacfunc =jac, guess = 17)
))

plot(sol2, type="l",lwd=2)
  
# including bound
bound <- function(i,y,parms){
  if (i ==1) return(y[1]-1)
  if (i ==2) return(y[2])
  if (i ==3) return(y[2])
}

print(system.time(
sol2b  <- bvptwp(bound = bound, leftbc = 2,x=x,
        func=mathieu2, jacfunc =jac, xguess = xguess,
        yguess = yguess)
))

### initial guess = 35 ####

xguess <-  c(0,1,2*pi)
yguess <- matrix(nr=3,rep(35,9))


print(system.time(
sol3  <- bvpshoot(bound = bound, leftbc = 2,x=x, atol=1e-9,
        func=mathieu2, jacfunc =jac, guess=c(y=1,dy=0,lambda=35))
))
        
## and jacbound
jacbound <- function(i,y,parms){
  if (i ==1) return(c(1,0,0))
  else return(c(0,1,0))
}
print(system.time(
sol3b  <- bvptwp(bound = bound, jacbound = jacbound, leftbc = 2, x=x,
        func=mathieu2, jacfunc =jac, xguess = xguess,
        yguess = yguess)
))

### initial guess = 105 ####

xguess <-  c(0,1,2*pi)
yguess <- matrix(nr=3,rep(105,9))

print(system.time(
sol4  <- bvpshoot(bound = bound, jacbound = jacbound, leftbc = 2,x=x,
        func=mathieu2, jacfunc =jac, guess=c(y=1,dy=1,lam=105))
))
print(system.time(
sol4b  <- bvptwp(bound = bound, jacbound = jacbound, leftbc = 2, x=x,
        func=mathieu2, jacfunc =jac, xguess = xguess,
        yguess = yguess)
))

par(mfrow=c(2,3))
plot(sol1,which="y", mfrow=NULL,type="l",lwd=2)
plot(sol2,which="y", mfrow=NULL,type="l",lwd=2)
plot(sol3,which="y", mfrow=NULL,type="l",lwd=2)
plot(sol1b,which="y", mfrow=NULL,type="l",lwd=2)
plot(sol2b,which=1, mfrow=NULL,type="l",lwd=2)
plot(sol3b,which=1, mfrow=NULL,type="l",lwd=2)
par(mfrow=c(1,1))


c(sol1[1,4],sol2[1,4],sol3[1,4],sol4[1,4])
c(sol1b[1,4],sol2b[1,4],sol3b[1,4],sol4b[1,4])
