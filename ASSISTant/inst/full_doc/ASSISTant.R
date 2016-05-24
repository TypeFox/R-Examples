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
data(LLL.SETTINGS)
str(LLL.SETTINGS)

## ------------------------------------------------------------------------
scenario <- LLL.SETTINGS$scenarios$S0
designParameters <- list(prevalence = LLL.SETTINGS$prevalences$table1,
                       mean = scenario$mean,
                       sd = scenario$sd)
designA <- ASSISTDesign$new(trialParameters = LLL.SETTINGS$trialParameters,
                            designParameters = designParameters)
print(designA)

## ------------------------------------------------------------------------
result <- designA$explore(numberOfSimulations = 5000, showProgress = FALSE)
analysis <- designA$analyze(result)
print(designA$summary(analysis))

## ------------------------------------------------------------------------
scenario <- LLL.SETTINGS$scenarios$S1
trueParameters <- list(prevalence = LLL.SETTINGS$prevalences$table1,
                       mean = scenario$mean,
                       sd = scenario$sd)
result <- designA$explore(numberOfSimulations = 5000,
                          trueParameters = trueParameters, showProgress = FALSE)
analysis <- designA$analyze(result)
print(designA$summary(analysis))

## ------------------------------------------------------------------------
scenario <- LLL.SETTINGS$scenarios$S2
trueParameters <- list(prevalence = LLL.SETTINGS$prevalences$table1,
                       mean = scenario$mean,
                       sd = scenario$sd)
result <- designA$explore(numberOfSimulations = 5000,
                          trueParameters = trueParameters, showProgress = FALSE)
analysis <- designA$analyze(result)
print(designA$summary(analysis))

## ------------------------------------------------------------------------
scenario <- LLL.SETTINGS$scenarios$S3
trueParameters <- list(prevalence = LLL.SETTINGS$prevalences$table1,
                       mean = scenario$mean,
                       sd = scenario$sd)
result <- designA$explore(numberOfSimulations = 5000,
                          trueParameters = trueParameters, showProgress = FALSE)
analysis <- designA$analyze(result)
print(designA$summary(analysis))

## ------------------------------------------------------------------------
scenario <- LLL.SETTINGS$scenarios$S4
trueParameters <- list(prevalence = LLL.SETTINGS$prevalences$table1,
                       mean = scenario$mean,
                       sd = scenario$sd)
result <- designA$explore(numberOfSimulations = 5000,
                          trueParameters = trueParameters, showProgress = FALSE)
analysis <- designA$analyze(result)
print(designA$summary(analysis))

## ------------------------------------------------------------------------
scenario <- LLL.SETTINGS$scenarios$S5
trueParameters <- list(prevalence = LLL.SETTINGS$prevalences$table1,
                       mean = scenario$mean,
                       sd = scenario$sd)
result <- designA$explore(numberOfSimulations = 5000,
                          trueParameters = trueParameters, showProgress = FALSE)
analysis <- designA$analyze(result)
print(designA$summary(analysis))

## ------------------------------------------------------------------------
scenario <- LLL.SETTINGS$scenarios$S6
trueParameters <- list(prevalence = LLL.SETTINGS$prevalences$table1,
                       mean = scenario$mean,
                       sd = scenario$sd)
result <- designA$explore(numberOfSimulations = 5000,
                          trueParameters = trueParameters, showProgress = FALSE)
analysis <- designA$analyze(result)
print(designA$summary(analysis))

## ------------------------------------------------------------------------
scenario <- LLL.SETTINGS$scenarios$S7
trueParameters <- list(prevalence = LLL.SETTINGS$prevalences$table1,
                       mean = scenario$mean,
                       sd = scenario$sd)
result <- designA$explore(numberOfSimulations = 5000,
                          trueParameters = trueParameters, showProgress = FALSE)
analysis <- designA$analyze(result)
print(designA$summary(analysis))

## ------------------------------------------------------------------------
scenario <- LLL.SETTINGS$scenarios$S8
trueParameters <- list(prevalence = LLL.SETTINGS$prevalences$table1,
                       mean = scenario$mean,
                       sd = scenario$sd)
result <- designA$explore(numberOfSimulations = 5000,
                          trueParameters = trueParameters, showProgress = FALSE)
analysis <- designA$analyze(result)
print(designA$summary(analysis))

## ------------------------------------------------------------------------
scenario <- LLL.SETTINGS$scenarios$S9
trueParameters <- list(prevalence = LLL.SETTINGS$prevalences$table1,
                       mean = scenario$mean,
                       sd = scenario$sd)
result <- designA$explore(numberOfSimulations = 5000,
                          trueParameters = trueParameters, showProgress = FALSE)
analysis <- designA$analyze(result)
print(designA$summary(analysis))

## ------------------------------------------------------------------------
scenario <- LLL.SETTINGS$scenarios$S10
trueParameters <- list(prevalence = LLL.SETTINGS$prevalences$table1,
                       mean = scenario$mean,
                       sd = scenario$sd)
result <- designA$explore(numberOfSimulations = 5000,
                          trueParameters = trueParameters, showProgress = FALSE)
analysis <- designA$analyze(result)
print(designA$summary(analysis))

## ------------------------------------------------------------------------
scenario <- LLL.SETTINGS$scenarios$S0
designParameters <- list(prevalence = LLL.SETTINGS$prevalences$table2,
                       mean = scenario$mean,
                       sd = scenario$sd)
designA <- ASSISTDesign$new(trialParameters = LLL.SETTINGS$trialParameters,
                            designParameters = designParameters)
print(designA)

## ------------------------------------------------------------------------
result <- designA$explore(numberOfSimulations = 5000, showProgress = FALSE)
analysis <- designA$analyze(result)
print(designA$summary(analysis))

## ------------------------------------------------------------------------
scenario <- LLL.SETTINGS$scenarios$S1
trueParameters <- list(prevalence = LLL.SETTINGS$prevalences$table2,
                       mean = scenario$mean,
                       sd = scenario$sd)
result <- designA$explore(numberOfSimulations = 5000, trueParameters = trueParameters,
                          showProgress = FALSE)
analysis <- designA$analyze(result)
print(designA$summary(analysis))

## ------------------------------------------------------------------------
scenario <- LLL.SETTINGS$scenarios$S2
trueParameters <- list(prevalence = LLL.SETTINGS$prevalences$table2,
                       mean = scenario$mean,
                       sd = scenario$sd)
result <- designA$explore(numberOfSimulations = 5000, trueParameters = trueParameters,
                          showProgress = FALSE)
analysis <- designA$analyze(result)
print(designA$summary(analysis))

## ------------------------------------------------------------------------
scenario <- LLL.SETTINGS$scenarios$S3
trueParameters <- list(prevalence = LLL.SETTINGS$prevalences$table2,
                       mean = scenario$mean,
                       sd = scenario$sd)
result <- designA$explore(numberOfSimulations = 5000, trueParameters = trueParameters,
                          showProgress = FALSE)
analysis <- designA$analyze(result)
print(designA$summary(analysis))

## ------------------------------------------------------------------------
scenario <- LLL.SETTINGS$scenarios$S4
trueParameters <- list(prevalence = LLL.SETTINGS$prevalences$table2,
                       mean = scenario$mean,
                       sd = scenario$sd)
result <- designA$explore(numberOfSimulations = 5000, trueParameters = trueParameters,
                          showProgress = FALSE)
analysis <- designA$analyze(result)
print(designA$summary(analysis))

## ------------------------------------------------------------------------
scenario <- LLL.SETTINGS$scenarios$S5
trueParameters <- list(prevalence = LLL.SETTINGS$prevalences$table2,
                       mean = scenario$mean,
                       sd = scenario$sd)
result <- designA$explore(numberOfSimulations = 5000, trueParameters = trueParameters,
                          showProgress = FALSE)
analysis <- designA$analyze(result)
print(designA$summary(analysis))

## ------------------------------------------------------------------------
scenario <- LLL.SETTINGS$scenarios$S6
trueParameters <- list(prevalence = LLL.SETTINGS$prevalences$table2,
                       mean = scenario$mean,
                       sd = scenario$sd)
result <- designA$explore(numberOfSimulations = 5000, trueParameters = trueParameters,
                          showProgress = FALSE)
analysis <- designA$analyze(result)
print(designA$summary(analysis))

## ------------------------------------------------------------------------
scenario <- LLL.SETTINGS$scenarios$S7
trueParameters <- list(prevalence = LLL.SETTINGS$prevalences$table2,
                       mean = scenario$mean,
                       sd = scenario$sd)
result <- designA$explore(numberOfSimulations = 5000, trueParameters = trueParameters,
                          showProgress = FALSE)
analysis <- designA$analyze(result)
print(designA$summary(analysis))

## ------------------------------------------------------------------------
scenario <- LLL.SETTINGS$scenarios$S8
trueParameters <- list(prevalence = LLL.SETTINGS$prevalences$table2,
                       mean = scenario$mean,
                       sd = scenario$sd)
result <- designA$explore(numberOfSimulations = 5000, trueParameters = trueParameters,
                          showProgress = FALSE)
analysis <- designA$analyze(result)
print(designA$summary(analysis))

## ------------------------------------------------------------------------
scenario <- LLL.SETTINGS$scenarios$S9
trueParameters <- list(prevalence = LLL.SETTINGS$prevalences$table2,
                       mean = scenario$mean,
                       sd = scenario$sd)
result <- designA$explore(numberOfSimulations = 5000, trueParameters = trueParameters,
                          showProgress = FALSE)
analysis <- designA$analyze(result)
print(designA$summary(analysis))

## ------------------------------------------------------------------------
scenario <- LLL.SETTINGS$scenarios$S10
trueParameters <- list(prevalence = LLL.SETTINGS$prevalences$table2,
                       mean = scenario$mean,
                       sd = scenario$sd)
result <- designA$explore(numberOfSimulations = 5000, trueParameters = trueParameters,
                          showProgress = FALSE)
analysis <- designA$analyze(result)
print(designA$summary(analysis))

