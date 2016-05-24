## duotrio <- function ()
## {
##   duotrio <- binomial()
##   duotrio$link <- "Link for the duo-trio test"
##   duotrio$linkinv <- function(eta) {
##     tres <- 1 - pnorm(eta/sqrt(2)) - pnorm(eta/sqrt(6)) +
##       2 * pnorm(eta/sqrt(2)) * pnorm(eta/sqrt(6))
##     tres[eta < 0] <- 0.5
##     tres
##   }
##   duotrio$mu.eta <- function(eta) -dnorm(eta/sqrt(2))/sqrt(2) -
##     dnorm(eta/sqrt(6))/sqrt(6) + 2 *
##       (dnorm(eta/sqrt(2)) * pnorm(eta/sqrt(6))/sqrt(2) +
##        pnorm(eta/sqrt(2)) * dnorm(eta/sqrt(6))/sqrt(6))
##   duotrio$linkfun <- function(mu) {
##     duotriog <- function(d, p) duotrio$linkinv(d) - p
##     tres <- mu
##     for(i in 1:length(mu)) {
##       if(mu[i] <= 0.5)
##         tres[i] <- 0
##       else if(mu[i] > 1-1e-9)
##         tres[i] <- Inf
##       else
##         tres[i] <- uniroot(duotriog, c(0, 14), p = mu[i])$root
##     }
##     tres
##   }
##   duotrio
## }

duotrio <- function() {
  duotrio <- binomial()
  duotrio$link <- "Link for the duo-trio test"
  duotrio$linkinv <- function(eta) {
    ok <- eta > 0 & eta < 20
    eta[eta <= 0] <- 0.5
    eta[eta >= 20] <- 1
    if(sum(ok)) {
      eta.ok <- eta[ok]
      pnorm.eta.2 <- pnorm(eta.ok * sqrt(1/2))
      pnorm.eta.6 <- pnorm(eta.ok * sqrt(1/6))
      eta[ok] <-
        1 - pnorm.eta.2 - pnorm.eta.6 + 2 * pnorm.eta.2 * pnorm.eta.6
    }
    pmin(pmax(eta, 0.5), 1) ## restrict to [0.5, 1] - just to be sure
  }
  duotrio$mu.eta <- function(eta) {
    ok <- eta > 0 ## no upper limit
    eta[!ok] <- 0
    if(sum(ok)) {
      eta.ok <- eta[ok]
      sqrt.2 <- sqrt(1/2)
      sqrt.6 <- sqrt(1/6)
      eta.2 <- eta.ok * sqrt.2
      eta.6 <- eta.ok * sqrt.6
      A <- dnorm(eta.2) * sqrt.2
      B <- dnorm(eta.6) * sqrt.6
      C <- dnorm(eta.2)
      D <- pnorm(eta.6) * sqrt.2
      E <- pnorm(eta.2)
      eta[ok] <- -A - B + 2 * (C * D + E * B)
    }
    pmax(eta, 0) ## gradient cannot be negative.
  }
  duotrio$linkfun <- function(mu) {
    eps <- 1e-10
    ok <- mu > 0.5 & mu < 1 - eps
    mu[mu <= 0.5] <- 0
    mu[mu >= 1 - eps] <- Inf
    if(sum(ok)) {
      duotriog <- function(d, p) duotrio$linkinv(d) - p
      mu[ok] <- sapply(mu[ok], function(mu) {
        uniroot(f=duotriog, interval=c(0, 16), p=mu)$root })
    }
    pmax(mu, 0) ## delta cannot be negative
  }
  duotrio
}

## `threeAFC` <-
##   function ()
## {
##   threeAFC <- binomial()
##   threeAFC$link <- "Link for the 3-AFC test"
##   threeAFC$linkinv <- function(eta) {
##     threeAFCg <- function(x, d) dnorm(x - d) * pnorm(x)^2
##     tres <- eta
##     for (i in 1:length(eta)) {
##       if (eta[i] > 0)
##         tres[i] <- integrate(threeAFCg, -Inf, Inf, d = eta[i])$value
##       else tres[i] <- 1/3
##     }
##     tres
##   }
##   threeAFC$mu.eta <- function(eta) {
##     threeAFCgd <- function(x, d) (x - d) * dnorm(x - d) *
##       pnorm(x)^2
##     tres <- eta
##     for (i in 1:length(eta)) tres[i] <-
##       integrate(threeAFCgd, -Inf, Inf, d = eta[i])$value
##     tres
##   }
##   threeAFC$linkfun <- function(mu) {
##     threeAFCg2 <- function(d, p) -p + threeAFC$linkinv(d)
##     tres <- mu
##     for (i in 1:length(mu)) {
##       if (mu[i] > 1/3)
##         res <- uniroot(threeAFCg2, c(0, 10), p = mu[i])$root
##       if (mu[i] <= 1/3)
##         res <- 0
##       tres[i] <- res
##     }
##     tres
##   }
##   threeAFC
## }

threeAFC <- function ()
{
  threeAFC <- binomial()
  threeAFC$link <- "Link for the 3-AFC test"
  threeAFC$linkinv <- function(eta) {
    ok <- eta > 0 & eta < 9
    eta[eta <= 0] <- 1/3
    eta[eta >= 9] <- 1
    if(sum(ok)) {
      threeAFCg <- function(x, d) dnorm(x - d) * pnorm(x)^2
      eta[ok] <- sapply(eta[ok], function(eta) {
        integrate(threeAFCg, -Inf, Inf, d = eta)$value })
    }
    pmin(pmax(eta, 1/3), 1) ## restrict to [1/3, 1] - just to be sure
  }
  threeAFC$mu.eta <- function(eta) {
    ok <- eta > 0 & eta < 9
    eta[eta <= 0] <- 0
    eta[eta >= 9] <- 0
### Note: integration is not reliable for eta > 9
    if(sum(ok)) {
      threeAFCgd <- function(x, d)
        (x - d) * dnorm(x - d) * pnorm(x)^2
      eta[ok] <- sapply(eta[ok], function(eta) {
        integrate(threeAFCgd, -Inf, Inf, d = eta)$value })
    }
    pmax(eta, 0) ## gradient cannot be negative.
  }
  threeAFC$linkfun <- function(mu) {
    eps <- 1e-8
    ok <- mu > 1/3 & mu < 1 - eps
    mu[mu <= 1/3] <- 0
    mu[mu >= 1 - eps] <- Inf
    if(sum(ok)) {
      threeAFCg2 <- function(d, p) threeAFC$linkinv(d) - p
      mu[ok] <- sapply(mu[ok], function(mu)
                       uniroot(threeAFCg2, c(0, 9), p = mu)$root)
    }
    pmax(mu, 0)
  }
  threeAFC
}

## `triangle` <-
##   function ()
## {
##   triangle <- binomial()
##     triangle$link <- "Link for the triangle test"
##   triangle$linkinv <- function(eta) {
##     triangleg <- function(x, d) 2 * dnorm(x) *
##       (pnorm(-x * sqrt(3) + d * sqrt(2/3)) +
##        pnorm(-x * sqrt(3) - d * sqrt(2/3)))
##     tres <- eta
##     for (i in 1:length(eta)) {
##       if (eta[i] > 0)
##         tres[i] <- integrate(triangleg, 0, Inf, d = eta[i])$value
##       else tres[i] <- 1/3
##     }
##     tres
##   }
##   triangle$mu.eta <- function(eta) {
##     trianglegd <- function(x, d) 2 * sqrt(2/3) * dnorm(x) *
##       (dnorm(-x * sqrt(3) + d * sqrt(2/3)) -
##        dnorm(-x * sqrt(3) - d * sqrt(2/3)))
##     tres <- eta
##     for (i in 1:length(eta)) tres[i] <-
##       integrate(trianglegd, 0, Inf, d = eta[i])$value
##     tres
##   }
##   triangle$linkfun <- function(mu) {
##     triangleg2 <- function(d, p) -p + triangle$linkinv(d)
##     tres <- mu
##     for (i in 1:length(mu)) {
##       if (mu[i] > 1/3)
##         tres[i] <- uniroot(triangleg2, c(0, 10), p = mu[i])$root
##       if (mu[i] <= 1/3)
##         tres[i] <- 0
##     }
##     tres
##   }
##   triangle
## }

triangle <- function()
{
  triangle <- binomial()
  triangle$link <- "Link for the triangle test"
  triangle$linkinv <- function(eta) {
    ok <- eta > 0 & eta < 20
    eta[eta <= 0] <- 1/3
    eta[eta >= 20] <- 1
    if(sum(ok))
      eta[ok] <-
        pf(q=3, df1=1, df2=1, ncp=eta[ok]^2*2/3, lower.tail=FALSE)
    pmin(pmax(eta, 1/3), 1) ## restrict to [1/3, 1] - just to be sure
  }
  triangle$mu.eta <- function(eta) {
    ok <- eta > 0 & eta < 20
    eta[eta <= 0] <- 0
    eta[eta >= 20] <- 0
    if(sum(ok)) {
      d <- eta[ok]
      eta[ok] <- sqrt(2/3) * dnorm(d/sqrt(6)) *
        (pnorm(d/sqrt(2)) - pnorm(-d/sqrt(2)))
    }
    pmax(eta, 0) ## gradient cannot be negative.
  }
  triangle$linkfun <- function(mu) {
    eps <- 1e-8
    ok <- mu > 1/3 & mu < 1 - eps
    mu[mu <= 1/3] <- 0
    mu[mu >= 1 - eps] <- Inf
    if(sum(ok)) {
      triangleg2 <- function(d, p) triangle$linkinv(d) - p
      mu[ok] <- sapply(mu[ok], function(mu)
                       uniroot(triangleg2, c(0, 15), p = mu)$root)
    }
    mu
  }
  triangle
}

tetrad <- function()
{
  tetrad <- binomial()
  tetrad$link <- "Link for the unspecified tetrad test"
  tetrad$linkinv <- function(eta) {
    eps <- 1e-8
    ok <- eta > eps & eta < 9
    eta[eta <= eps] <- 1/3
    eta[eta >= 9] <- 1
    ## TetradsLinkinv <- function(eta) {
    ## ok <- eta > -9 & eta < 9
    ## eta[eta <= 0] <- 1/3
    ## eta[eta >= 9] <- 1
    if(sum(ok)) {
      tetrads.fun <- function(z, delta)
        dnorm(z) * (2 * pnorm(z) * pnorm(z - delta) -
                    pnorm(z - delta)^2)
      eta[ok] <- sapply(eta[ok], function(eta) {
        1 - 2*integrate(tetrads.fun, -Inf, Inf, delta=eta)$value })
    }
    pmin(pmax(eta, 1/3), 1) ## restrict to [1/3, 1] - just to be sure
  }
  tetrad$mu.eta <- function(eta) {
    eps <- 1e-8
    ok <- eta > eps & eta < 9
    eta[eta <= eps] <- 0
    eta[eta >= 9] <- 0
     if(sum(ok)) {
      Linkinv <- function(eta) {
        tetrads.fun <- function(z, delta)
          dnorm(z) * (2 * pnorm(z) * pnorm(z - delta) -
                      pnorm(z - delta)^2)
        sapply(eta, function(eta) {
          1 - 2*integrate(tetrads.fun, -Inf, Inf, delta=eta)$value })
      }
      eta[ok] <- sapply(eta[ok], function(eta) grad(Linkinv, eta))
### FIXME: Could probably do the integration by hand here.
    }
    pmax(eta, 0) ## gradient cannot be negative.
  }
  tetrad$linkfun <- function(mu) {
    eps <- 1e-8 ## What is the right eps here?
    ok <- mu > 1/3 & mu < 1 - eps
    mu[mu <= 1/3] <- 0
    mu[mu >= 1 - eps] <- Inf
    if(sum(ok)) {
      tetrads <- function(d, p) tetrad$linkinv(d) - p
      mu[ok] <- sapply(mu[ok], function(mu)
                       uniroot(tetrads, c(0, 9), p = mu)$root)
    }
    pmax(mu, 0)
  }
  tetrad
}

## twoAFC <- function() {
##   twoAFC <- binomial()
##   twoAFC$link <- "Link for the 2-AFC test"
##   twoAFC$linkfun <- function(mu) {
##     tres <- mu
##     for (i in 1:length(mu)) {
##       if (mu[i] > 0.5)
##         tres[i] <- sqrt(2) * qnorm(mu[i])
##       if (mu[i] <= 0.5)
##         tres[i] <- 0
##     }
##     tres
##   }
##   twoAFC$linkinv <- function(eta) {
##     tres <- pnorm(eta/sqrt(2))
##     tres[eta < 0] <- 0.5
##     tres
##   }
##   twoAFC$mu.eta <- function(eta) dnorm(eta/sqrt(2))/sqrt(2)
##   twoAFC
## }

twoAFC <- function() {
  twoAFC <- binomial()
  twoAFC$link <- "Link for the 2-AFC test"
  twoAFC$linkinv <- function(eta) {
    ok <- eta > 0
    eta[!ok] <- 0.5
    eta[ok] <- pnorm(eta[ok] / sqrt(2))
    pmin(pmax(eta, 0.5), 1) ## restrict to [0.5, 1] - just to be sure
  }
  twoAFC$mu.eta <- function(eta) {
    ok <- eta >= 0
    eta[!ok] <- 0
    if(any(ok)) {
      sqrt.2 <- sqrt(1/2)
      eta[ok] <- dnorm(eta[ok] * sqrt.2) * sqrt.2
    }
    pmax(eta, 0) ## gradient cannot be negative.
  }
  twoAFC$linkfun <- function(mu) {
    ok <- mu > 0.5 & mu < 1
    mu[mu <= 0.5] <- 0
    mu[mu >= 1] <- Inf
    mu[ok] <- sqrt(2) * qnorm(mu[ok])
    pmax(mu, 0) ## delta cannot be negative
  }
  twoAFC
}
