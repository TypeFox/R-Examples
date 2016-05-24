test.convergence.density <- function()
{
  s0 <- 100
  d0 <- .2
  mu <- 0.2
  sigmaS <- 0.3
  kappa <- 0.9
  alpha <- 0.03
  sigmaE <- 0.3
  rho <- -0.5
  lambda <- 0.1
  time = 3
  Time = 5.01
  r = 0.2


  obj <- schwartz2f(s0 = s0, delta0 = d0, mu = mu, sigmaS = sigmaS,
                    alpha = alpha, kappa = kappa, sigmaE = sigmaE, rho =rho)
  r.fut <- rfutures(1000000, time = time, ttm = Time, obj, lambda = lambda, r = r, measure = "Q")

  g0 <- pricefutures(Time, obj, r = r,lambda = lambda)
  K <- g0 * 1.3
  ## For European call options
  p.opt <- priceoption(time = time, Time = Time, K = K, g0 = g0, sigmaS = sigmaS, kappa = kappa,
                       sigmaE = sigmaE, rho = rho, r = r, type = "call")
  payoff <- r.fut - K
  payoff[payoff < 0] <- 0
  p.opt.sim <- mean(payoff) * exp(-r * time)

  checkTrue(abs((p.opt - p.opt.sim) / p.opt) < 0.01, "Call: Monte carlo options valuation must converge")

  ## For European put options
  p.opt <- priceoption(time = time, Time = Time, K = K, g0 = g0, sigmaS = sigmaS, kappa = kappa,
                       sigmaE = sigmaE, rho = rho, r = r, type = "put")
  payoff <- K - r.fut
  payoff[payoff < 0] <- 0
  p.opt.sim <- mean(payoff) * exp(-r * time)

  checkTrue(abs((p.opt - p.opt.sim) / p.opt) < 0.01, "Put: Monte carlo options valuation must converge")

}

test.convergence.monte.carlo <- function()
{

  s0 <- 88
  d0 <- .2
  mu <- 0.2
  sigmaS <- 0.15
  kappa <- 1.9
  alpha <- 0.0343
  sigmaE <- 0.443
  rho <- 0.55
  lambda <- 0.18
  time = 3.11
  Time = 4.44
  r = 0.223


  obj <- schwartz2f(s0 = s0, delta0 = d0, mu = mu, sigmaS = sigmaS,
                    alpha = alpha, kappa = kappa, sigmaE = sigmaE, rho =rho)
  g0 <- pricefutures(Time, obj, r = r,lambda = lambda)
  K <- g0 * 1.11


  ## For European call options
  int.fun <- function(x, obj, t, T, K, lambda, r){
    (x - K) * dfutures(x, time = t, ttm = T, obj, lambda = lambda, r = r, measure = "Q")
  }
  p.opt.num <- integrate(int.fun, lower = K, upper = Inf, obj = obj, t = time, T = Time, K = K,
                         lambda = lambda, r = r)$value * exp(-r * time)

  p.opt <- priceoption(time = time, Time = Time, K = K, g0 = g0, sigmaS = sigmaS, kappa = kappa,
                       sigmaE = sigmaE, rho = rho, r = r, type = "call")

  checkTrue(abs((p.opt - p.opt.num) / p.opt) < 0.01,
            "Call: Integration method must give same result as analytical solution")


  ## For European put options
  int.fun <- function(x, obj, t, T, K, lambda, r){
    (K - x) * dfutures(x, time = t, ttm = T, obj, lambda = lambda, r = r, measure = "Q")
  }
  p.opt.num <- integrate(int.fun, lower = 0, upper = K, obj = obj, t = time, T = Time, K = K,
                         lambda = lambda, r = r)$value * exp(-r * time)

  p.opt <- priceoption(time = time, Time = Time, K = K, g0 = g0, sigmaS = sigmaS, kappa = kappa,
                       sigmaE = sigmaE, rho = rho, r = r, type = "put")

  checkTrue(abs((p.opt - p.opt.num) / p.opt) < 0.01,
            "Put: Integration method must give same result as analytical solution")
}
