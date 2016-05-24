## Copyright (C) 2013 Marius Hofert, Bernhard Pfaff
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


momest <- function(data, trials, limit = 10){
  out <- rep(NA, limit)
  for(k in 1:limit){
    term <- 1.0
    for(j in 1:k){
      term <- (term * (data - (j - 1.0)))/(trials - (j - 1.0))
    }
    out[k] <- mean(term)
  }
  out
}
##
cal.probitnorm <-function(pi1, pi2){
  rootfunc <- function(x, k1, k2){
    k2 - pmvnorm(lower=-Inf, upper=c(qnorm(k1), qnorm(k1)), corr = equicorr(2, x))
  }
  out <- uniroot(rootfunc, c(0.0001, 0.999), k1 = pi1, k2 = pi2)
  rho <- out$root
  mu <- qnorm(pi1) / sqrt(1 - rho)
  sigma <- sqrt(rho) / sqrt(1 - rho)
  c(mu = mu, sigma = sigma, rho.asset = rho)
}
##
cal.beta <- function(pi1, pi2){
  a <- (pi1 * (pi1 - pi2))/(pi2 - pi1^2.)
  b <- ((1. - pi1) * (pi1 - pi2))/(pi2 - pi1^2.)
  c(a = a, b = b)
}
##
cal.claytonmix <- function(pi1, pi2){
  rootfunc <- function(theta, k1, k2){
    log(2 * k1^( - theta) - 1) + theta * log(k2)
  }
  out <- uniroot(rootfunc, c(1e-010, 1), k1 = pi1, k2 = pi2)
  theta <- out$root
  c(pi = pi1, theta = theta)
}
##
pprobitnorm <- function(q, mu, sigma){
  pnorm((qnorm(q) - mu)/sigma)
}
##
pclaytonmix <- function(q, pi, theta){
  1 - pgamma((( - log(q))/(pi^( - theta) - 1)), (1/theta))
}
##
dprobitnorm <- function(x, mu, sigma){
  (dnorm((qnorm(x) - mu) / sigma)) / (dnorm(qnorm(x)) * sigma)
}
##
dclaytonmix <- function(x, pi, theta){
  dgamma((( - log(x))/(pi^( - theta) - 1)), (1/theta))/(x * (pi^( - theta) - 1))
}
##
rprobitnorm <- function(n, mu, sigma){
  pnorm(rnorm(n, mu, sigma))
}
##
rlogitnorm <- function(n, mu, sigma){
  (1 + exp(-rnorm(n, mu, sigma)))^(-1)
}
##
rclaytonmix <- function(n, pi, theta){
  exp(-(pi^(-theta) - 1) * rgamma(n, 1 / theta))
}
##
rtcopulamix <- function(n, pi, rho.asset, df){
  W <- df / rchisq(n, df)
  THETA <- rnorm(n)
  pnorm((qt(pi, df)/sqrt(W) - THETA * sqrt(rho.asset)) / sqrt(1 - rho.asset))
}
##
rbinomial.mixture <- function(n = 1000, m = 100, model = c("probitnorm", "logitnorm", "beta"), ...){
  model <- match.arg(model)
  mixdist <- eval(parse(text = paste("r", model, sep = "")))
  Q <- mixdist(n, ...)
  rbinom(n, m, Q)
}
##
fit.binomial <- function(M, m){
  phat <- sum(M)/sum(m)
  loglik <- function(p, M, m){
    logb(p) * sum(M) + logb(1. - p) * (sum(m) - sum(M))
  }
  llmax <- loglik(phat, M, m)
  pse <- sqrt((phat * (1 - phat)) / length(M))
  pi <- phat
  pi2 <- phat^2
  rhoY <- 0
  list(par.ests = phat, par.ses = pse, maxloglik = llmax, pi = pi, pi2 = pi2, rhoY = rhoY)
}
##
fit.binomialBeta <- function(M, m, startvals = c(2, 2), ses = FALSE, ...){
  negloglik <- function(theta, defaults, trials){
    length(trials) * base::lbeta(theta[1]^2, theta[2]^2) - sum(base::lbeta(theta[1]^2 + defaults, b = theta[2]^2 + trials - defaults))
  }
  truenegloglik <- function(theta, defaults, trials){
    length(trials) * base::lbeta(theta[1], theta[2]) - sum(base::lbeta(theta[1] + defaults, theta[2] + trials - defaults))
  }
  theta <- sqrt(startvals)
  fit <- optim(theta, negloglik, defaults = M, trials = m, ...)
  par.ests <- fit$par
  ifelse(fit$convergence == 0, conv <- TRUE, conv <- FALSE)
  loglik <-  - negloglik(par.ests, defaults=M, trials=m)
  par.ests <- par.ests^2
  if(ses){
    hessmatrix <- hessian(truenegloglik, par.ests, defaults=M, trials=m)
    varcov <- solve(hessmatrix)
    par.ses <- sqrt(diag(varcov))
  } else {
    par.ses <- rep(NA, 2)
  }
  pi <- par.ests[1]/sum(par.ests)
  pi2 <- (par.ests[1] + 1) / (sum(par.ests) + 1) * pi
  rhoY <- (pi2 - pi^2) / (pi - pi^2)
  list(par.ests = par.ests, par.ses = par.ses, maxloglik = loglik, converged = conv, pi = pi, pi2 = pi2, rhoY = rhoY, fit = fit)
}
##
fit.binomialProbitnorm <- function(M, m, startvals = c(-1, 0.5), ...){
  link <- pnorm
  integrand <- function(p, kk, n, mu, sigma, lfunc){
    exp(kk * log(lfunc(mu + sigma * qnorm(p))) + (n - kk) * log(1 - lfunc(mu + sigma * qnorm(p))))
  }
  integrand2 <- function(p, mu, sigma, kk, lfunc){
    exp(kk * log(lfunc(mu + sigma * qnorm(p))))
  }
  MLdata <- cbind(M,m)
  negloglik <- function(theta, defaultData, FNIntegrate = integrand, FNlink = link){
    n <- dim(defaultData)[1]
    out <- numeric(n)
    for(i in 1:n){
      out[i] <-  - logb(integrate(FNIntegrate, lower = 0.0, upper = 1.0,  kk = defaultData[i, 1], n = defaultData[i, 2], mu = theta[1], sigma = theta[2]^2, lfunc = FNlink)$value)
    }
    sum(out)
  }
  theta <- c(startvals[1], sqrt(startvals[2]))
  fit <- nlminb(theta, negloglik, defaultData = MLdata, FNIntegrate = integrand, FNlink = link, ...)
  par.ests <- fit$par
  ifelse(fit$convergence == 0, conv <- TRUE, conv <- FALSE)
  loglik <-  -fit$objective
  par.ests[2] <- par.ests[2]^2
  tmp1 <- integrate(f = integrand2, lower = 0, upper = 1, mu = par.ests[1], sigma = par.ests[2], kk = 1, lfunc = link)
  pi <- tmp1$value
  tmp2 <- integrate(f = integrand2, lower = 0, upper = 1,  mu = par.ests[1], sigma = par.ests[2], kk = 2, lfunc = link)
  pi2 <- tmp2$value
  rhoY <- (pi2 - pi^2)/(pi - pi^2)
  list(par.ests = par.ests, maxloglik = loglik, converged = conv, pi = pi, pi2 = pi2, rhoY = rhoY, fit = fit)
}
##
fit.binomialLogitnorm <- function(M, m, startvals = c(-1, 0.5), ...){
  link <- function(z){
    1/(1 + exp( - z))
  }
  integrand <- function(p, kk, n, mu, sigma, lfunc){
    exp(kk * log(lfunc(mu + sigma * qnorm(p))) + (n - kk) * log(1 - lfunc(mu + sigma * qnorm(p))))
  }
  integrand2 <- function(p, mu, sigma, kk, lfunc){
    exp(kk * log(lfunc(mu + sigma * qnorm(p))))
  }
  MLdata <- cbind(M,m)
  negloglik <- function(theta, defaultData, FNIntegrate = integrand, FNlink = link){
    n <- dim(defaultData)[1]
    out <- numeric(n)
    for(i in 1:n){
      out[i] <-  - logb(integrate(FNIntegrate, lower = 0.0, upper = 1.0,  kk = defaultData[i, 1], n = defaultData[i, 2], mu = theta[1], sigma = theta[2]^2, lfunc = FNlink)$value)
    }
    sum(out)
  }
  theta <- c(startvals[1], sqrt(startvals[2]))
  fit <- nlminb(theta, negloglik, defaultData = MLdata, FNIntegrate = integrand, FNlink = link, ...)
  par.ests <- fit$par
  ifelse(fit$convergence == 0, conv <- TRUE, conv <- FALSE)
  loglik <-  -fit$objective
  par.ests[2] <- par.ests[2]^2
  tmp1 <- integrate(f = integrand2, lower = 0, upper = 1, mu = par.ests[1], sigma = par.ests[2], kk = 1, lfunc = link)
  pi <- tmp1$value
  tmp2 <- integrate(f = integrand2, lower = 0, upper = 1,  mu = par.ests[1], sigma = par.ests[2], kk = 2, lfunc = link)
  pi2 <- tmp2$value
  rhoY <- (pi2 - pi^2)/(pi - pi^2)
  list(par.ests = par.ests, maxloglik = loglik, converged = conv, pi = pi, pi2 = pi2, rhoY = rhoY, fit = fit)
}
