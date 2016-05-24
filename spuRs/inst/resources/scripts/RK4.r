# Code spuRs/resources/scripts/RK4.r

RK4 <- function(dydt, t0, y0, h, n, ...) {
  # Runge-Kutta 4 applied to the system of ODEs with grad ftn dydt
  # and initial condition y0 at time t0
  # simulates n steps of size h, returned in a matrix with n + 1 rows 
  # and cols corresponding to the elements of y
  # additional arguments are passed to dydt
  y <- matrix(0, nrow = n+1, ncol = length(y0))
  y[1,] <- y0
  t <- t0
  for (i in 1:n) {
    k1 <- dydt(t, y[i,], ...)
    k2 <- dydt(t + h/2, y[i,] + h*k1/2, ...)
    k3 <- dydt(t + h/2, y[i,] + h*k2/2, ...)
    k4 <- dydt(t + h, y[i,] + h*k3, ...)
    y[i+1,] <- y[i,] + h*(k1/6 + k2/3 + k3/3 + k4/6)
    t <- t + h
  }
  return(y)
}

