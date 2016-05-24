# Code spuRs/resources/scripts/midpoint.r

midpoint <- function(dydt, t0, y0, h, n, ...) {
  # midpoint scheme applied to the system of ODEs with grad ftn dydt
  # and initial condition y0 at time t0
  # simulates n steps of size h, returned in a matrix with n + 1 rows 
  # and cols corresponding to the elements of y
  # additional arguments are passed to dydt
  y <- matrix(0, nrow = n+1, ncol = length(y0))
  y[1,] <- y0
  t <- t0
  for (i in 1:n) {
    y2 <- y[i,] + h/2*dydt(t, y[i,], ...)
    y[i+1,] <- y[i,] + h*dydt(t + h/2, y2, ...)
    t <- t + h
  }
  return(y)
}

