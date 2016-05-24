# Code spuRs/resources/scripts/RK4adapt.r
# loadable spuRs function

RK4adapt <- function(dydt, t0, y0, t1, h0 = 1, tol = 1e-10, ...) {
  # 4-th order Runge-Kutta with adaptive step size
  # dydt(t, y, ...) gives the gradient of y(t)
  # t0 and y0 are starting conditions
  # The system is solved up to time t1. Initial step size is h0, 
  # which is adjusted so that the error in each step is approx tol.
  # Returns a list with elements t a vector giving times,
  # and y a matrix whose rows give the soltn at successive times.  
  t <- t0  # current time
  y <- y0  # current state
  h <- h0  # current step size
  t.vec <- t
  y.vec <- matrix(y, nrow=1)
  while (t < t1) {
    # RK4 using step size h
    k1 <- dydt(t, y, ...)
    k2 <- dydt(t + h/2, y + h*k1/2, ...)
    k3 <- dydt(t + h/2, y + h*k2/2, ...)
    k4 <- dydt(t + h, y + h*k3, ...)
    y1 <- y + h*(k1/6 + k2/3 + k3/3 + k4/6)
    # RK4 using two steps size h/2
    k1a <- k1
    k2a <- dydt(t + h/4, y + h/2*k1a/2, ...)
    k3a <- dydt(t + h/4, y + h/2*k2a/2, ...)
    k4a <- dydt(t + h/2, y + h/2*k3a, ...)
    y2 <- y + h/2*(k1a/6 + k2a/3 + k3a/3 + k4a/6)
    k1b <- dydt(t + h/2, y2, ...)
    k2b <- dydt(t + 3*h/4, y2 + h/2*k1b/2, ...)
    k3b <- dydt(t + 3*h/4, y2 + h/2*k2b/2, ...)
    k4b <- dydt(t + h, y2 + h/2*k3b, ...)
    y3 <- y2 + h/2*(k1b/6 + k2b/3 + k3b/3 + k4b/6)
    # update h (decrease)
    if (max(abs(y3 - y1)) > tol) {
      h <- h/2
    } else {
      # update t and y
      t <- t + h
      y <- y3 + (y3 - y1)/15
      # update h (increase)
      if (max(abs(y3 - y1)) < tol/2) h <- 3*h/2
      # record keeping
      t.vec <- c(t.vec, t)
      y.vec <- rbind(y.vec, y)
    }
  }
  return(list(t=t.vec, y=y.vec))
}


