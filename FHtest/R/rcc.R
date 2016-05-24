rcc <-
function (times, status, z, rho, lambda)
{
  uz <- sort(unique(z))
  k <- length(uz)
  fit <- survfit(Surv(times, status) ~ 1)

  fail <- fit$time[fit$n.event > 0]

  w <- fit$surv[fit$n.event > 0]
  w <- c(1, w[1:(length(w) - 1)])

  neventg <- table(z[status > 0], times[status > 0])
  nevent <- colSums(neventg)

  nriskg <- matrix(1, length(fail), k)
  for (i in 1:k)
    nriskg[, i] <- colSums(matrix(rep(fail,each=sum(z==uz[i])),,length(fail))<=times[z==uz[i]])
  nrisk <- rowSums(nriskg)

  observed <- (w^rho * (1 - w)^lambda) %*% t(neventg)
  expected <- (w^rho * (1 - w)^lambda) %*% (nriskg * (nevent/nrisk))

  v <- (w^rho * (1 - w)^lambda)^2 * (nevent * (nrisk - nevent)/(nrisk - 1))
  v[nrisk == 1] <- 0
  v <- diag(c(v %*% (nriskg/nrisk))) - t(nriskg/nrisk) %*% ((nriskg/nrisk) * v)

  list(obs = c(observed), exp = c(expected), var = v)
}
