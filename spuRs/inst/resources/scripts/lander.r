# Code spuRs/resources/scripts/lander.r

lander <- function(Ti, a, h = 0.1, tol = 1e-2, plot.traj = FALSE) {
  # Solve the trajectory of a rocket.
  # The rocket falls until time Ti (s), then produces thrust a (kg m/s^2), 
  # until the velocity is +ve, it hits the ground, or the fuel runs out, 
  # at which point the simulation stops.
  #
  # h is initial step size for adaptive RK4
  # tol is tolerance for adaptive RK4
  # plot.traj flag for plotting trajectories
  #
  # state y = (velocity m/s, height m, mass kg)
  # lander returns the final value of y
  
  v0 = -1000   # initial velocity m/s
  ht0 = 50000  # initial height m
  m0 <- 10000  # unladen mass of rocket kg
  f <- 1500    # initial fuel load kg
  r <- 0.0001  # mass of fuel burnt per unit thrust per unit time s/m
  g <- 1.633   # acceleration due to gravity m/s^2
  
  # derivatives of y at time t
  dydt <- function(t, y) {
    if (t < Ti) {
      return(c(-g, y[1], 0))
    } else {
      return(c(a/y[3] - g, y[1], -r*a))
    }
  }
  
  # 4-th order Runge-Kutta with adaptive step size
  t <- 0
  y <- c(v0, ht0, m0 + f)
  if (plot.traj) y.trace <- c(t, y)
  while (y[1] < 0 && y[2] > 0 && y[3] > m0) { 
    # readline(paste(t, y[1], y[2]))  # for step-by-step reporting
    # RK4 using step size h
    k1 <- dydt(t, y)
    k2 <- dydt(t + h/2, y + h*k1/2)
    k3 <- dydt(t + h/2, y + h*k2/2)
    k4 <- dydt(t + h, y + h*k3)
    y1 <- y + h*(k1/6 + k2/3 + k3/3 + k4/6)
    # RK4 using two steps size h/2
    k1a <- k1
    k2a <- dydt(t + h/4, y + h/2*k1a/2)
    k3a <- dydt(t + h/4, y + h/2*k2a/2)
    k4a <- dydt(t + h/2, y + h/2*k3a)
    y2 <- y + h/2*(k1a/6 + k2a/3 + k3a/3 + k4a/6)
    k1b <- dydt(t + h/2, y2)
    k2b <- dydt(t + 3*h/4, y2 + h/2*k1b/2)
    k3b <- dydt(t + 3*h/4, y2 + h/2*k2b/2)
    k4b <- dydt(t + h, y2 + h/2*k3b)
    y3 <- y2 + h/2*(k1b/6 + k2b/3 + k3b/3 + k4b/6)
    y.new <- y3 + (y3 - y1)/15
    # update h (decrease)
    if (max(abs(y3 - y1)) > tol || 
        y.new[1] >= tol || y.new[2] <= -tol || y.new[3] <= m0 - tol) {
      h <- h/2
    } else {
      # update t and y
      t <- t + h
      y <- y.new
      # update h (increase)
      if (max(abs(y3 - y1)) < tol/2) h <- 3*h/2
      # record keeping
      if (plot.traj) y.trace <- rbind(y.trace, c(t, y))
    }
  }  

  # output
  if (plot.traj) {
    opar <- par(mfrow=c(3,1), mar=c(4,4,1,1))
    plot(y.trace[,1], y.trace[,2], type='o', xlab='t', ylab='velocity')
    plot(y.trace[,1], y.trace[,3], type='o', xlab='t', ylab='height')
    plot(y.trace[,1], y.trace[,4], type='o', xlab='t', ylab='mass')
    par <- opar
  }
  return(y)
}
