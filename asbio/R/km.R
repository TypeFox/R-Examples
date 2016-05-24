km <- function (r, d, var = "O", conf = 0.95, age.seq = seq(1, length(r)), 
                ylab = "Survivorship", xlab = "Age class", 
                ...) 
{
  step <- (1 - (d/r))
  s.hat <- matrix(nrow = length(r), ncol = 1)
  s.hat[1] <- step[1]
  for (i in 2:length(r)) {
    s.hat[i] <- s.hat[i-1] * step[i]
  }
  
  Green.Sum <- matrix(nrow = length(r), ncol = 1)
  Green.Var <- matrix(nrow = length(r), ncol = 1)
  for (i in 1:length(r)) {
    Green.Var[1] <- 0
    Green.Sum[1] <- d[1]/(r[1] * (r[1] - d[1]))
    Green.Sum[i + 1] <- d[i + 1]/(r[i + 1] * (r[i + 1] - 
                                                d[i + 1])) + Green.Sum[i]
    Green.Var[i + 1] <- s.hat[i + 1]^2 * Green.Sum[i + 1]
    Green.Var <- Green.Var[-(length(r) + 1)]
  }
  Oakes.Var <- (s.hat^2 * (1 - s.hat)/r)
  if (var == "O") {
    C.L <- s.hat - qnorm(1 - ((1 - conf)/2)) * ((Oakes.Var)^0.5)
    C.U <- s.hat + qnorm(1 - ((1 - conf)/2)) * ((Oakes.Var)^0.5)
  }
  if (var == "G") {
    C.L <- s.hat - qnorm(1 - ((1 - conf)/2)) * ((Green.Var)^0.5)
    C.U <- s.hat + qnorm(1 - ((1 - conf)/2)) * ((Green.Var)^0.5)
  }
  plot(age.seq, s.hat, type = "b", xlab = xlab, ylab = ylab, 
       ...)
  lines(age.seq, C.L, lty = 2)
  lines(age.seq, C.U, lty = 2)
  res <- list()
  res$s.hat <- s.hat
  res$Greenwood.Var <- Green.Var
  res$Oakes.Var <- Oakes.Var
  res$CI <- cbind(C.L, C.U)
  res
}
