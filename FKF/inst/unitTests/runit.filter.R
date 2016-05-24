test.compare.to.stats <- function()
{

  ## Fit a l

  data <- as.numeric(treering)
  data[c(3, 10:31)] <- NA

  fit.stats <- StructTS(data, type = "level")

  yt <- rbind(data)
  a0 <- yt[1]
  P0 <- matrix(100)
  dt <- ct <- matrix(0)
  Zt <- Tt <- matrix(1)
  HHt <- var(t(yt), na.rm = TRUE)
  GGt <- var(t(yt), na.rm = TRUE) * 0.1

  StructTS.level <- function(par, ...){
    return(-fkf(HHt = matrix(par[1]), GGt = matrix(par[2]), ...)$logLik)
  }
  fit.fkf <- optim(c(HHt, GGt), StructTS.level, yt = yt,
                   a0 = a0, P0 = P0, dt = dt, ct = ct,
                   Zt = Zt, Tt = Tt, check.input = FALSE)

  checkTrue(all((fit.fkf$par / fit.stats$coef - 1) < 0.01),
            "Difference between 'FKF' and 'stats' implementation is greater than 1 percent!")
}
