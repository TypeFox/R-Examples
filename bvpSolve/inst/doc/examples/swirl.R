## =============================================================================
## "Swirling Flow III", a test problem in Ascher, Mattheij,
## and Russell, Numerical Solution of Boundary Value Problems
## for Ordinary Differential Equations", Classics in Applied
## Mathematics Series, SIAM, Philadelphia, 1995].
## g'' = (g f' - f g')/eps
## f'''' = (-ff'''-gg')/eps
## g(0)=-1,f(0)=0,f'(0)=0, g(1)=1,f(1)=0,f'(1)=0
##
## 1st order system (y1=g, y3=f)
## y1' = y2
## y2' = (y1*y4 -y3*y2)/eps
## y3'=y4
## y4'=y5
## y5'=y6
## y6'=(-y3y6-y1y2)/eps
## y1(0)=-1,y3(0)=0,y4(0)=0, y1(1)=1,y3(1)=0,y4(1)=0
## =============================================================================

require(bvpSolve)

# Declare global problem dependent parameters.
x       <- seq(0, 1, 0.01)
yini <- c(-1, NA, 0, 0, NA, NA)
yend <- c(1 , NA, 0, 0, NA, NA)


fsub <- function (t, Y, eps)  {
 return(list(c(f1 <- Y[2],
                f2 <- (Y[1]*Y[4] - Y[3]*Y[2])/eps,
	              f3 <- Y[4],
              	f4 <- Y[5],
              	f5 <- Y[6],
	              f6 <- (-Y[3]*Y[6] - Y[1]*Y[2])/eps)))
}

# Solve the model. Simple call
# shooting does not work with eps this small.

# bvptwp does...
print(system.time(Sol <- bvptwp(x = x, func = fsub, 
                  yini = yini, yend = yend, eps = 0.01, parms = 0.01)))

pairs(Sol, col = "blue", lty = 2)

# For successively smaller values of eps:
print(system.time(Sol2 <- bvptwp(x = x, func = fsub, 
                  xguess = Sol[,1], yguess = t(Sol[,-1]),
                  yini = yini, yend = yend, eps=0.001, epsini = 0.01,
                  parms=0.001)))
pairs(Sol2, col = "red", lty = 2)

# For successively smaller values of eps:
print(system.time(Sol3 <- bvptwp(x = x, func = fsub, 
                  xguess = Sol2[,1], yguess = t(Sol2[,-1]),
                  yini = yini, yend = yend, eps = 0.0001, 
                  epsini = 0.001, parms = 0.0001)))
pairs(Sol3, col = "darkgreen", lty = 2)


# For successively smaller values of eps:
print(system.time(Sol4 <- bvptwp(x = x, func = fsub, 
                  xguess = Sol3[,1], yguess = t(Sol3[,-1]),
                  yini = yini, yend = yend, eps = 5e-5, 
                  epsini = 1e-4, parms = 5e-5)))
pairs(Sol4, col = "orange", lty = 2)

# with conditioning 
print(system.time(Sol5 <- bvptwp(atol = 1e-5, x = x, func = fsub, cond = TRUE, 
                  xguess = Sol3[,1], yguess = t(Sol3[,-1]),
                  yini = yini, yend = yend, eps = 5e-5 , parms=5e-5)))
pairs(Sol5, col = "yellow2", lty = 2)
