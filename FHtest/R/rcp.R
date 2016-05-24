rcp <-
function (times, ind, rho, lambda)
{
  sv <- survfit(Surv(times, ind) ~ 1)
  w <- c(1, sv$surv[1:(length(sv$surv) - 1)])
  w <- w^rho * (1 - w)^lambda
  aux <- -w * sv$n.event/sv$n.risk
  aux <- cumsum(aux)
  aux.event <- (aux + w)[sv$n.event > 0]
  aux.censor <- aux[sv$n.censor > 0]
  cc <- numeric(length(ind))
  cc[ind == 1] <- aux.event[match(times[ind == 1], sv$time[sv$n.event > 0])]
  cc[ind == 0] <- aux.censor[match(times[ind == 0], sv$time[sv$n.censor > 0])]
  cc
}
