# Code spuRs/resources/scripts/euler.r

euler <- function(dydt, t0, y0, h, n, ...) {
  # Euler's method applied to the system of ODEs with grad ftn dydt
  # and initial condition y0 at time t0
  # simulates n steps of size h, returned in a matrix with n + 1 rows 
  # and cols corresponding to the elements of y
  # additional arguments are passed to dydt
  y <- matrix(0, nrow = n + 1, ncol = length(y0))
  y[1,] <- y0
  t <- t0
  for (i in 1:n) {
    y[i+1,] <- y[i,] + h*dydt(t, y[i,], ...)
    t <- t + h
  }
  return(y)
}
