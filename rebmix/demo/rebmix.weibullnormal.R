##############################################
## R sources for reproducing the results in ##
##   Mitja Franko, Marko Nagode:            ##
##   Probability Density Function of the    ##
##   Equivalent Stress Amplitude using      ## 
##   Statistical Transformation             ##
##############################################

options(prompt = "> ", continue = "+ ", width = 70,
  useFancyQuotes = FALSE, digits = 3)

###########################
## Weibullnormal dataset ##
###########################

data("weibullnormal", package = "rebmix")

########## FlexMix ##########

library("flexmix")

set.seed(8)

# Weibull-normal mixture model.

WNmodel <- function (formula = .~.) {    
  retval <- new("FLXMC", weighted = TRUE,
    formula = formula, dist = "Weibull-normal",
    name = "Weibull-normal parameters")

  retval@defineComponent <- expression({
    logLik <- function(x, y) {
      logLik <- array(0)

      for (i in 1:nrow(y)) {
        logLik[i] <- dweibull(y[i, 1], shape = beta, scale = theta, log = TRUE) + 
                     dnorm(y[i, 2], mean = mean, sd = sd, log = TRUE)
      }

      logLik
    }

    new("FLXcomponent",
      parameters = list(mean = mean, sd = sd, beta = beta, theta = theta), df = df, logLik = logLik)
  })

  retval@fit <- function(x, y, w, ...) {
    n <- nrow(y)

    mean <- sum(w * y[, 2]) / sum(w)

    sd <- sqrt(sum(w * (y[, 2] - mean)^2) / sum(w))
      
    beta <- 0.9

    i <- 1; Error <- 1
    while((i <= 1000) && (Error == 1)) {
      T0 <- log(y[, 1])
      T1 <- exp(T0 * beta)

      A0 <- sum(w * T0)
      A1 <- sum(w * T1 * T0)
      A2 <- sum(w * T1)
      A3 <- sum(w * T1 * T0 * T0)
      A4 <- sum(w)

      f <- (1.0 / beta) + A0 / A4 - A1 / A2
      df <- ((A1 / A2) * (A1 / A2)) - (A3 / A2) - (1.0 / beta / beta)
        
      dbeta <- f / df

      beta <- beta - dbeta
        
      if (is.nan(dbeta) || is.infinite(dbeta) || (beta <= 1e-6)) {
        stop("Not converged!")
      }

      if (abs(dbeta / beta) < 1e-6) {
        Error <- 0
      }

      i <- i + 1
    }

    theta <- exp(log(A2 / A4) * (1.0 / beta))

    df <- 2 * ncol(y)

    para <- list(sd = sd, mean = mean, beta = beta, theta = theta, n.obs = n)
 
    with(para, eval(retval@defineComponent))
  }
 
  retval
} ## WNmodel

# Estimate number of components, component weights and component parameters by FlexMix.

timeFlexMix <- system.time(weibullnormalestFlexMix <- flexmix(as.matrix(weibullnormal) ~ 1, k = 3, model = WNmodel()))

parameters(weibullnormalestFlexMix)
prior(weibullnormalestFlexMix)
timeFlexMix
summary(weibullnormalestFlexMix)

########## REBMIX ##########

library("rebmix")

# Estimate number of components, component weights and component parameters.

n = nrow(weibullnormal)

Sturges <- as.integer(1 + log2(n)) # Minimum v follows the Sturges rule.
Log10 <- as.integer(10 * log10(n)) # Maximum v follows the Log10 rule.

timerebmix <- system.time(weibullnormalestrebmix <- REBMIX(Dataset = 
  list(weibullnormal = weibullnormal),
  Preprocessing = "histogram", cmax = 5, Criterion = "AIC",
  pdf = c("Weibull", "normal"),
  K = kseq(Sturges, Log10, 0.07)))

weibullnormalestrebmix
AIC(weibullnormalestrebmix)
timerebmix 

plot(weibullnormalestrebmix, nrow = 2, ncol = 3, 
  what = c("density", "marginal", "IC", "logL"), npts = 1000)
  
summary(weibullnormalestrebmix)

detach("package:rebmix")
detach("package:flexmix")