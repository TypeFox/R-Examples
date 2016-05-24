## ----echo=F--------------------------------------------------------------
### get knitr just the way we like it

knitr::opts_chunk$set(
  message = FALSE,
  warning = FALSE,
  error = FALSE,
  tidy = FALSE,
  cache = FALSE
)

## ------------------------------------------------------------------------
library(ASSISTant)
##Fix randomization vector N, errors, eps
trialParameters <- list(N = c(200, 340, 476), type1Error = 0.025,
                        eps = 1/2, type2Error = 0.1)

## ------------------------------------------------------------------------
designParameters <- list(
    nul0 = list(prevalence = rep(1/6, 6), mean = matrix(0, 2, 6),
                sd = matrix(1, 2, 6)),
    alt1 = list(prevalence = rep(1/6, 6), mean = rbind(rep(0, 6),
                                                       c(0.5, 0.4, 0.3, 0, 0, 0)),
                sd = matrix(1, 2, 6)),
    alt2 = list(prevalence = rep(1/6, 6), mean = rbind(rep(0, 6),
                                                     c(0.5, 0.5, 0, 0, 0, 0)),
                sd = matrix(1,2, 6)),
    alt3 = list(prevalence = rep(1/6, 6), mean = rbind(rep(0, 6), rep(0.36, 6)),
                sd = matrix(1,2, 6)),
    alt4 = list(prevalence = rep(1/6, 6), mean = rbind(rep(0, 6), rep(0.30, 6)),
                sd = matrix(1,2, 6)),
    alt5 = list(prevalence = rep(1/6, 6), mean = rbind(rep(0, 6),
                                                       c(0.4, 0.3, 0.2, 0, 0, 0)),
                sd = matrix(1,2, 6)),
    alt6 = list(prevalence = rep(1/6, 6), mean = rbind(rep(0, 6),
                                                       c(0.5, 0.5, 0.3, 0.3, 0.1, 0.1)),
                sd = matrix(1,2, 6))
)

## ------------------------------------------------------------------------
defuse3 <- DEFUSE3Design$new(trialParameters = trialParameters,
                             designParameters = designParameters$nul0, showProgress = FALSE)
print(defuse3)

## ------------------------------------------------------------------------
result <- defuse3$explore(numberOfSimulations = 5000, showProgress = FALSE)
analysis <- defuse3$analyze(result)
print(defuse3$summary(analysis))

## ------------------------------------------------------------------------
result <- defuse3$explore(numberOfSimulations = 5000, trueParameters = designParameters$alt1,
                          showProgress = FALSE)
analysis <- defuse3$analyze(result)
print(defuse3$summary(analysis))

## ------------------------------------------------------------------------
result <- defuse3$explore(numberOfSimulations = 5000, trueParameters = designParameters$alt2,
                          showProgress = FALSE)
analysis <- defuse3$analyze(result)
print(defuse3$summary(analysis))

## ------------------------------------------------------------------------
result <- defuse3$explore(numberOfSimulations = 5000, trueParameters = designParameters$alt3,
                          showProgress = FALSE)
analysis <- defuse3$analyze(result)
print(defuse3$summary(analysis))

## ------------------------------------------------------------------------
result <- defuse3$explore(numberOfSimulations = 5000, trueParameters = designParameters$alt4,
                          showProgress = FALSE)
analysis <- defuse3$analyze(result)
print(defuse3$summary(analysis))

## ------------------------------------------------------------------------
result <- defuse3$explore(numberOfSimulations = 5000, trueParameters = designParameters$alt5,
                          showProgress = FALSE)
analysis <- defuse3$analyze(result)
print(defuse3$summary(analysis))

## ------------------------------------------------------------------------
result <- defuse3$explore(numberOfSimulations = 5000, trueParameters = designParameters$alt6,
                          showProgress = FALSE)
analysis <- defuse3$analyze(result)
print(defuse3$summary(analysis))

