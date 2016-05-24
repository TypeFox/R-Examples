chute <- function(Ti, h = 0.1, tol = 1e-4, plot.traj = F) {
  # Missing code is denoted ??
  #
  # Trajectory of a skydiver
  # Ti is time to open parachute
  # state y = (velocity, height)

  g <- 9.81   # acceleration due to gravity m/s^2
  m <- 100    # mass of skydiver kg
  d1 <- 0.31  # drag without chute kg/m
  d2 <- 56    # drag with chute kg/m
  ht0 <- 2000 # jump height m
  # d1 is obtained from terminal velocity for a 100kg person of 200 km/h
  # d2 is obtained from terminal velocity for a 100kg parachute of 15 km/h
  
  # derivatives of y at time x
  dydx <- function(x, y) {
    ??
  }
  
  # 4-th order Runge-Kutta with adaptive step size
  x <- 0
  y <- ??  # initial conditions
  if (plot.traj) y.trace <- c(x, y)
  while (??) {  # stopping criteria
    # readline(paste(x, y[1], y[2]))  # for diagnostics
    # RK4 using step size h
    k1 <- h*dydx(x, y)
    k2 <- h*dydx(x + h/2, y + k1/2)
    k3 <- h*dydx(x + h/2, y + k2/2)
    k4 <- h*dydx(x + h, y + k3)
    y1 <- y + k1/6 + k2/3 + k3/3 + k4/6
    # RK4 using two steps size h/2
    k1a <- k1/2
    k2a <- h/2*dydx(x + h/4, y + k1a/2)
    k3a <- h/2*dydx(x + h/4, y + k2a/2)
    k4a <- h/2*dydx(x + h/2, y + k3a)
    y2 <- y + k1a/6 + k2a/3 + k3a/3 + k4a/6
    k1b <- h/2*dydx(x + h/2, y2)
    k2b <- h/2*dydx(x + 3*h/4, y2 + k1b/2)
    k3b <- h/2*dydx(x + 3*h/4, y2 + k2b/2)
    k4b <- h/2*dydx(x + h, y2 + k3b)
    y3 <- y2 + k1b/6 + k2b/3 + k3b/3 + k4b/6
    y.new <- y3 + (y3 - y1)/15
    # update h
    if (max(abs(y3 - y1)) > tol || y.new[2] < -tol) {
      h <- h/2
    } else {
      # update x and y
      x <- x + h
      y <- y.new
      # update h
      if (max(abs(y3 - y1)) < tol/2) h <- 3*h/2
      # record keeping
      if (plot.traj) y.trace <- rbind(y.trace, c(x, y))
    }
  }  
  
  # output
  if (plot.traj) {
    opar <- par(mfrow=c(2,1), mar=c(4,4,1,1))
    plot(y.trace[,1], y.trace[,2], type='o', xlab='t', ylab='velocity')
    plot(y.trace[,1], y.trace[,3], type='o', xlab='t', ylab='height')
    par <- opar
  }
  return(c(x, y))
}

# safe landing speed is 20 km/h = 5.6 m/s
