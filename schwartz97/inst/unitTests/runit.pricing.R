test.futures.relations <- function()
  {

    object <- schwartz2f(s0 = 1.1, delta0 = 0.2, alpha = 0.1)
    checkEquals(pricefutures(ttm = 0, object, r = 0, lambda = 0), 1.1,
                "Futures must converge to spot price")

    checkTrue(pricefutures(ttm = 0.5, object, r = 0, lambda = 0.1) < 1.1,
              "Positive market price of convenience yield risk must lead to lower futures price")

  }

test.futures.valuation <- function()
  {

    lambda <- -0.2
    r <- 0.14
    obj <- schwartz2f(rho = -0.3, sigmaS = .8, sigmaE = 0.5)

    ## Time to maturity = 1.3
    ttm <- 1.3
    r.fut <- pricefutures(ttm, obj, lambda = lambda, r = r)
    r.fut.sample <- mean(rfutures(1e6, ttm, ttm, obj, lambda = lambda, r = r, measure = "Q"))
    checkTrue(abs((r.fut - r.fut.sample) / r.fut) < 0.01,
              "ttm = 1.3: Simulated futures prices must converge to analytical price!")

    ## Time to maturity = 3.3
    ttm <- 3.3
    r.fut <- pricefutures(ttm, obj, lambda = lambda, r = r)
    r.fut.sample <- mean(rfutures(1e6, ttm, ttm, obj, lambda = lambda, r = r, measure = "Q"))
    checkTrue(abs((r.fut - r.fut.sample) / r.fut) < 0.01,
              "ttm = 3.3: Simulated futures prices must converge to analytical price!")


    ## Time to maturity = 6.9
    ttm <- 6.9
    r.fut <- pricefutures(ttm, obj, lambda = lambda, r = r)
    r.fut.sample <- mean(rfutures(5e6, ttm, ttm, obj, lambda = lambda, r = r, measure = "Q"))
    checkTrue(abs((r.fut - r.fut.sample) / r.fut) < 0.03,
              "ttm = 6.9: Simulated futures prices must converge to analytical price!")

  }
