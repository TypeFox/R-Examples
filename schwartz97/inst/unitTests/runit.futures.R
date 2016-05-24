test.convergence <- function()
{
  object <- schwartz2f(s0 = 4.4, alpha = 0.4, delta0 = -0.2, rho = -.8)
  r.state <- rstate(100000, time = 2.2, object, method = "chol")[,1]
  r.fut <- rfutures(100000, time = 2.2, ttm = 2.2, object, r = 0)
  checkTrue(abs(mean(r.state) - mean(r.fut)) < 0.05,
            "Futures must converge to spot")
}
